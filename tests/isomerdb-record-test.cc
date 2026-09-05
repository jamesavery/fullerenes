// Correctness tests for IsomerDB::Entry::from_dual (src/c++/isomerdb_record.cc),
// the record producer that replaced the Fortran program's.
//
// The isomer database is the oracle: every record in it was written by the
// legacy program from the same isomer, so rebuilding the dual from a record's
// own spiral and producing the record again must give the record back, field
// for field (HLgap to the column's rounding).  The error codes are exercised
// on inputs that provoke them where such an input can be constructed:
// NotRepresentable from a cage with more faces than an 8-bit column holds.
// NoRegularSpiral needs a fullerene with no jump-free spiral at all (none
// below C380), and UnknownPointGroup a graph Symmetry cannot classify as one
// of the 28 fullerene groups -- neither has a small constructible witness.

#include "fullerene-test-main.hh"

#include "fullerenes/isomerdb.hh"
#include "fullerenes/triangulation.hh"

#include <cstdio>
#include <vector>

using Entry = IsomerDB::Entry;
using Field = IsomerDB::Field;
using EntryResult = IsomerDB::EntryResult;

namespace {

const std::vector<int> dodecahedron_rspi{0,1,2,3,4,5,6,7,8,9,10,11};

}  // namespace

// Every record of every requested size, All corpus: the dual rebuilt from the
// record's spiral produces the record.
TEST(IsomerRecord, ReproducesEveryDatabaseRecord) {
  ASSERT_FALSE(fullerene_test::sizes().empty()) << "--sizes= would make this sweep assert nothing";
  size_t n_compared = 0;
  for (int N : fullerene_test::sizes()) {
    if (!IsomerDB::is_installed(N)) continue;
    IsomerDB db = IsomerDB::readPDB(N, false);
    ASSERT_GT(db.entries.size(), 0u) << "empty database for C" << N;
    const Field stored = db.stored_fields();
    for (size_t i = 0; i < db.entries.size(); i++) {
      const Entry& e = db.entries[i];
      const FullereneDual D(N, e.rspi_zero_based());
      const EntryResult r = Entry::from_dual(D, stored);
      ASSERT_TRUE(r.ok()) << "C" << N << " isomer " << i + 1 << ": " << r.name() << ": " << r.why;
      const Field differ = r.entry.diff(e, stored, Entry::HLgap_column_tol);
      EXPECT_EQ(differ, Field::None) << "C" << N << " isomer " << i + 1 << " differs in "
                                     << IsomerDB::field_names(differ) << "\n  db: " << e.to_string()
                                     << "\n  produced: " << r.entry.to_string();
      n_compared++;
    }
    printf("[          ] C%d: %zu records reproduced\n", N, db.entries.size());
  }
  if (n_compared == 0) GTEST_SKIP() << "none of the requested sizes is installed";
}

// The producer mask is honoured: fields outside it stay zero, the key is
// always produced.
TEST(IsomerRecord, ProducesOnlyTheRequestedFields) {
  const FullereneDual D(20, dodecahedron_rspi);
  const EntryResult key_only = Entry::from_dual(D, Field::None);
  ASSERT_TRUE(key_only.ok()) << key_only.why;
  EXPECT_EQ(key_only.entry.rspi_zero_based(), dodecahedron_rspi);
  EXPECT_EQ(key_only.entry.NeHOMO, 0);
  EXPECT_EQ(key_only.entry.ncycham, 0);
  EXPECT_EQ(std::string(key_only.entry.group, 3), std::string(3, '\0'));

  const EntryResult all = Entry::from_dual(D, Field::All);
  ASSERT_TRUE(all.ok()) << all.why;
  EXPECT_EQ(all.entry.point_group().to_string(), "Ih");
  EXPECT_EQ(all.entry.ncycham, 30);
  EXPECT_EQ(all.entry.NeHOMO, 2);         // 18 electrons fill 3, sqrt5 x3, 1 x5; two in the
  EXPECT_EQ(all.entry.NedgeHOMO, 4);      // 4-fold zero level: an open shell
  EXPECT_EQ(all.entry.diff(key_only.entry, Field::RSPI), Field::None);
}

// A cage with more faces than the RSPI column's 255: the record cannot exist,
// and the producer says so instead of storing a truncated position.
TEST(IsomerRecord, TooManyFacesIsNotRepresentable) {
  const FullereneDual small(20, dodecahedron_rspi);
  const FullereneDual big(small.GCtransform(5, 1));   // C620, 312 faces
  ASSERT_EQ(big.N, 312);
  const EntryResult r = Entry::from_dual(big, Field::None);
  EXPECT_EQ(r.code, EntryResult::Code::NotRepresentable) << r.name() << ": " << r.why;
  EXPECT_FALSE(r.why.empty());
}

int main(int argc, char** argv) {
  return fullerene_test::run(argc, argv, {20, 24, 26, 28, 30, 32, 34, 36, 38, 40}, "isomerdb-record");
}
