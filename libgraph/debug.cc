#include "debug.hh"

map<string,int> Debug::channel_verbosity;
const char*     Debug::verbosity_level_txt[5] = {"ERROR","WARNING","INFO","INFO","INFO"};

ostream Debug::nullstream(0);
int     Debug::prefix = 1;

map<string,int>      MathematicaDebug::channel_verbosity;
map<string,ostream*> MathematicaDebug::channel_stream;

ostream MathematicaDebug::nullstream(0);
int     MathematicaDebug::prefix = 0;
