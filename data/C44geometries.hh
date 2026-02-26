#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

// Binary file layout:
//   c044all.rspi.u8:       12 uint8_t per record (0-indexed pentagon positions)
//   c044all.geometry.f32:  44*3 float32 per record (optimized fullerene coordinates)
//
// Record count = file_size / record_size.

struct C44Geometries {
  int n_isomers;
  std::vector<uint8_t>  rspi;   // n_isomers * 12
  std::vector<float>    geom;   // n_isomers * 44 * 3

  const uint8_t* pentagon_indices(int i) const { return &rspi[i * 12]; }
  const float*   points(int i)          const { return &geom[i * 44 * 3]; }

  static C44Geometries load(const char* rspi_path, const char* geom_path) {
    C44Geometries db;

    FILE* fr = fopen(rspi_path, "rb");
    FILE* fg = fopen(geom_path, "rb");
    if(!fr || !fg){
      fprintf(stderr, "C44Geometries::load: cannot open %s or %s\n", rspi_path, geom_path);
      abort();
    }

    // Determine record count from rspi file size
    fseek(fr, 0, SEEK_END);
    long rspi_size = ftell(fr);
    fseek(fr, 0, SEEK_SET);
    db.n_isomers = (int)(rspi_size / 12);

    db.rspi.resize(db.n_isomers * 12);
    db.geom.resize(db.n_isomers * 44 * 3);

    if(fread(db.rspi.data(), 1, db.rspi.size(), fr) != db.rspi.size() ||
       fread(db.geom.data(), sizeof(float), db.geom.size(), fg) != db.geom.size()){
      fprintf(stderr, "C44Geometries::load: short read\n");
      abort();
    }

    fclose(fr);
    fclose(fg);
    return db;
  }
};
