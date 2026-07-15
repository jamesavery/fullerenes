// IsomerDB robustness tests: the readers must throw informative errors on a
// missing or corrupt database (never abort or return a poisoned sentinel),
// and read correctly when the database is present.
//
// The failure cases redirect IsomerDB::database_path to a nonexistent root
// and restore it afterwards, so they are independent of the configured
// FULLERENE_DATABASE_PATH; the happy-path cases use the configured database
// and skip (rather than fail) when it is absent, so the suite also runs on
// machines without the database installed.

#include <gtest/gtest.h>
#include <fstream>
#include "fullerenes/isomerdb.hh"

namespace {

// RAII redirect of the static database root.
struct DatabasePathGuard {
  std::string saved;
  explicit DatabasePathGuard(const std::string& path) : saved(IsomerDB::database_path) {
    IsomerDB::database_path = path;
  }
  ~DatabasePathGuard() { IsomerDB::database_path = saved; }
};

bool database_present() {
  std::ifstream f(IsomerDB::database_path + "/All/c060all.database");
  return f.good();
}

}  // namespace

TEST(IsomerDBErrors, ReadPDBMissingDatabaseThrows) {
  DatabasePathGuard guard("/nonexistent-fullerene-database");
  try {
    IsomerDB::readPDB(60, false);
    FAIL() << "readPDB on a missing database must throw";
  } catch (const std::runtime_error& e) {
    // The message must name the file and how to fix the configuration.
    EXPECT_NE(std::string(e.what()).find("c060all.database"), std::string::npos) << e.what();
    EXPECT_NE(std::string(e.what()).find("FULLERENE_DATABASE_PATH"), std::string::npos) << e.what();
  }
}

TEST(IsomerDBErrors, ReadBinaryMissingDatabaseThrows) {
  DatabasePathGuard guard("/nonexistent-fullerene-database");
  EXPECT_THROW(IsomerDB::readBinary(60, false), std::runtime_error);
}

TEST(IsomerDBErrors, GetIsomerMissingDatabaseThrows) {
  DatabasePathGuard guard("/nonexistent-fullerene-database");
  EXPECT_THROW(IsomerDB::getIsomer(60, 1, false), std::runtime_error);
}

TEST(IsomerDBErrors, NumberIsomersOutOfRangeIsZero) {
  EXPECT_EQ(IsomerDB::number_isomers(18), 0);      // below C20
  EXPECT_EQ(IsomerDB::number_isomers(21), 0);      // odd N
  EXPECT_EQ(IsomerDB::number_isomers(58, "Any", true), 0);   // below IPR base C60
  EXPECT_EQ(IsomerDB::number_isomers(60), 1812);   // table lookup, no file I/O
  EXPECT_EQ(IsomerDB::number_isomers(60, "Any", true), 1);
}

TEST(IsomerDBReads, ReadPDB_C60) {
  if (!database_present()) GTEST_SKIP() << "fullerene database not installed";
  IsomerDB db = IsomerDB::readPDB(60, false);
  EXPECT_EQ(db.N, 60);
  EXPECT_EQ(db.Nisomers, 1812);
  EXPECT_EQ((int)db.entries.size(), 1812);
}

TEST(IsomerDBReads, MakeIsomerFromEntry) {
  if (!database_present()) GTEST_SKIP() << "fullerene database not installed";
  IsomerDB db = IsomerDB::readPDB(60, false);
  FullereneGraph g = IsomerDB::makeIsomer(60, db.entries[0]);
  EXPECT_EQ(g.N, 60);
}
