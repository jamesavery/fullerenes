#ifndef DEBUG_HH
#define DEBUG_HH

#include "auxiliary.hh"

class Debug {
public:
  static map<string,int> channel_verbosity; 
  enum {ERROR,WARNING,INFO1,INFO2,INFO3} verbosity_level;
  static const char* verbosity_level_txt[5];

  string channel;
  int level, prefix;

  static ostream nullstream;

  Debug(string channel, int level, int prefix=1) : channel(channel), level(level), prefix(prefix) {}

  template <typename T> friend ostream& operator<<(const Debug& S, const T& x)
  {
    if(channel_verbosity[S.channel] < S.level && channel_verbosity["all"] < S.level) return nullstream;
    if(S.prefix>1) cerr << S.channel << "::";
    if(S.prefix>0) cerr << verbosity_level_txt[S.level] << ": ";
    cerr << x;

    return cerr;
  }
};

class MathDebug {
public:
  static map<string,int>       channel_verbosity; 
  static map<string,ostream*>  channel_stream;

  string channel;
  int level;
  static int prefix;

  static ostream nullstream;
 
  MathDebug(string channel, int level, int prefix=1) : channel(channel), level(level) {}

  template <typename T> friend ostream& operator<<(const MathDebug& S, const T& x)
  {
    if(channel_verbosity[S.channel] < S.level && channel_verbosity["all"] < S.level) return nullstream;

    ostream& stream = *channel_stream[S.channel];

    stream << x;
    stream.flush();

    return stream;
  }
};

#endif
