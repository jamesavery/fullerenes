#include <sqlite3.h>
#include "libgraph/isomerdb.hh"

string bigint(const int *digits, int base, int length)
{
  for(int i=0;
}

bool sqlite_execute_void_once(sqlite3 *database, string query_string)
{
  sqlite3_stmt *query = 0;
  int ret = 0;
  if((ret = sqlite3_prepare_v2(database,query_string.c_str(),query_string.size(),&query,0))!=SQLITE_OK){
    cerr << "Couldn't parse query: '"<<query_string<<"': " << sqlite3_errmsg(database) << "\n";
    return false;
  }
  if(sqlite3_step(query) != SQLITE_DONE){
    cerr << "Couldn't execute query: '"<<query_string<<"': " << sqlite3_errmsg(database) << " \n";
    return false;
  }
  sqlite3_finalize(query);

  return true;
}

bool sqlite_insert_entry(sqlite3 *database, int id, IsomerDB::Entry e)
{
  ostringstream valuestring("");
  valuestring << "(" 
	      << id << ",'"
	      << string(e.group,3) << "','"
	      <<vector<int>(e.RSPI,e.RSPI+12)<<"','"
	      <<vector<int>(e.PNI,e.PNI+5)<<"','"
	      <<vector<int>(e.HNI,e.HNI+5)<<"',"
	      <<e.NeHOMO<<","<<e.NedgeHOMO<<","<<e.HLgap << ","
	      <<e.ncycham<<",'"
	      <<vector<int>(e.INMR,e.INMR+6)<<"')"
    ;

  return sqlite_execute_void_once(database,"insert into ENTRIES "
				  "(ID,SYMGROUP,RSPI,PNI,HNI,NEHOMO,NEDGEHOMO,HLGAP,NCYCHAM,INMR) "
				  "values " + valuestring.str());
}

int main(int ac, char **av)
{
  if(ac<=2) return -1;
  
  int N         = strtol(av[1],0,0);
  string outdir = string(av[2]);

  IsomerDB DB(IsomerDB::readPDB(N));

  string filename = outdir+"/c"+pad_string(to_string(N),3)+(DB.IPR?"IPR":"all")+".sqlite";

  sqlite3 *database = 0;

  sqlite3_initialize();
  if(sqlite3_open(filename.c_str(),&database) != SQLITE_OK){
    cerr << "Couldn't create database on " << filename << ".\n";
    return -2;
  }

  string createtable = "create table ENTRIES ("
    "ID integer primary key, "
    "SYMGROUP  char(3), "
    "RSPI      tinyint(12), "
    "PNI       tinyint(5), "
    "HNI       tinyint(6), "
    "NEHOMO    tinyint, "
    "NEDGEHOMO tinyint, "
    "HLGAP     float, "
    "NCYCHAM   int, "
    "INMR      tinyint(6))";

  string createindex = "create unique index RSPIINDEX on ENTRIES ( RSPI )";

  if(sqlite_execute_void_once(database,"drop table if exists entries") == false) return -3;
  if(sqlite_execute_void_once(database,createtable) == false) return -3;
  //  if(sqlite_execute_void_once(database,createindex) == false) return -4;

  for(int i=0;i<DB.Nisomers;i++)
    sqlite_insert_entry(database,i+1,DB.entries[i]);

  return 0;
}
