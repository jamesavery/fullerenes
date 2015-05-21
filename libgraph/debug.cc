#include "debug.hh"

map<string,int> Debug::channel_verbosity;
const char* Debug::verbosity_level_txt[5] = {"ERROR","WARNING","INFO","INFO","INFO"};
ostream Debug::nullstream(0);
