#include "fullerenes/isomerdb.hh"
#include "fullerenes/symmetry.hh"
#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/hueckel.hh"

// IsomerDB::Entry::from_dual -- the database record of an isomer, produced
// from its graph by the C++ producers that replace the Fortran ones:
//
//   RSPI                    FullereneDualView::regular_rspi (the database key)
//   PNI, HNI                FullereneDualView::neighbour_indices (spiral.f DualAnalyze)
//   group, NMR              Symmetry::point_group / NMR_pattern
//   NeHOMO, NedgeHOMO,      hueckel::analyze (hueckel.f HueckelAnalyze +
//   HLgap                     the gap-zeroing convention of spiral.f's writer)
//   ncycham                 GraphView::hamilton_cycle_count (hamilton.f HamiltonCyc)
//
// Every way the record can fail to exist is a named EntryResult code returned
// at the point it is discovered; representability is checked here, before the
// with_* words whose contracts assume it, so nothing is thrown or caught.

using Entry = IsomerDB::Entry;
using Field = IsomerDB::Field;
using EntryResult = IsomerDB::EntryResult;

namespace {

using Code = EntryResult::Code;

bool fits_u8(long v) { return 0 <= v && v <= 255; }

// A value a record column cannot hold makes the record NotRepresentable, not
// truncated; the message names the column and the value.
EntryResult not_representable(const std::string& column, long v, const char* width = "8-bit") {
  return EntryResult::error(Code::NotRepresentable,
                            column + " = " + std::to_string(v) + " does not fit its " + width + " record column");
}

}  // namespace

EntryResult Entry::from_dual(const FullereneDualView& D, Field producers)
{
  Entry e{};
  const std::vector<int> rspi = D.regular_rspi();
  if (rspi.empty())
    return EntryResult::error(Code::NoRegularSpiral,
                              "the isomer has no jump-free ring spiral, so no database key");
  for (int p : rspi) if (!fits_u8(p + 1)) return not_representable("a face position", p + 1);
  e = with_rspi(e, rspi);

  if (!!(producers & (Field::PNI | Field::HNI)))
    e = with_neighbour_indices(e, D.neighbour_indices());

  if (!!(producers & Field::FromPointGroup)) {
    Symmetry S(D);
    const PointGroup pg = S.point_group();
    if (!pg.is_fullerene_group())
      return EntryResult::error(Code::UnknownPointGroup,
                                "Symmetry classified the isomer as '" + pg.to_string() +
                                "', not one of the 28 fullerene point groups");
    e = with_point_group(e, pg);
    // (count, orbit size) ascending in size -- the record's INMR order.
    const std::vector<std::pair<int,int>> nmr = S.NMR_pattern();
    if (nmr.size() > 3)
      return EntryResult::error(Code::NotRepresentable,
                                "the NMR pattern has " + std::to_string(nmr.size()) +
                                " orbit sizes; the record holds 3");
    for (size_t k = 0; k < nmr.size(); k++) {
      const auto [count, size] = nmr[k];
      if (!fits_u8(count)) return not_representable("an NMR orbit count", count);
      if (!fits_u8(size))  return not_representable("an NMR orbit size", size);
      e.INMR[2*k]   = u_int8_t(count);
      e.INMR[2*k+1] = u_int8_t(size);
    }
  }

  if (!!(producers & (Field::FromSpectrum | Field::Ncycham))) {
    const FullereneGraph G(D.dual_graph());   // the cage itself: 3-regular, stride 3
    if (!!(producers & Field::FromSpectrum)) {
      const hueckel::Analysis H = hueckel::analyze(G);
      if (!fits_u8(H.homo_electrons()))  return not_representable("the HOMO electron count", H.homo_electrons());
      if (!fits_u8(H.homo_degeneracy())) return not_representable("the HOMO degeneracy", H.homo_degeneracy());
      e.NeHOMO    = u_int8_t(H.homo_electrons());
      e.NedgeHOMO = u_int8_t(H.homo_degeneracy());
      e.HLgap     = float(H.gap);
    }
    if (!!(producers & Field::Ncycham)) {
      const int64_t h = G.hamilton_cycle_count();
      if (h < 1 || h > Entry::ncycham_max) return not_representable("the Hamilton-cycle count", h, "I7");
      e.ncycham = int(h);
    }
  }
  return {Code::Ok, e, ""};
}
