#ifndef DEBUG_HH
#define DEBUG_HH

#include "auxiliary.hh"

class Debug {
public:
  static map<string,int> channel_verbosity;

  string channel;
  int level;

  enum {ERROR,WARNING,INFO1,INFO2,INFO3} verbosity_level;
  static const char* verbosity_level_txt[5];

  Debug(string channel, int level) : channel(channel), level(level) {}

  template <typename T> friend ostream& operator<<(const Debug& S, const T& x)
  {
    if(channel_verbosity[S.channel] < S.level) return cerr;
    cerr << S.channel << "::" << verbosity_level_txt[S.level] << ": " << x;
    return cerr;
  }
};

#endif
