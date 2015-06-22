#ifndef DEBUG_HH
#define DEBUG_HH

#include "auxiliary.hh"

class Debug {
public:
  static map<string,int> channel_verbosity; 
  enum {ERROR,WARNING,INFO1,INFO2,INFO3} verbosity_level;
  static const char* verbosity_level_txt[5];

  string channel;
  int level;

  static ostream nullstream;
  static int prefix;

  Debug(string channel, int level) : channel(channel), level(level) {}

  template <typename T> friend ostream& operator<<(const Debug& S, const T& x)
  {
    if(channel_verbosity[S.channel] < S.level && channel_verbosity["all"] < S.level) return nullstream;
    if(prefix>1) cerr << S.channel << "::";
    if(prefix>0) cerr << verbosity_level_txt[S.level] << ": ";
    cerr << x;

    return cerr;
  }
};

class MathematicaDebug {
public:
  static map<string,int>      channel_verbosity; 
  static map<string,ostream*> channel_stream;

  string channel;
  int level;

  static ostream nullstream;
  static int prefix;

  MathematicaDebug(string channel, int level=1) : channel(channel), level(level) {}

  template <typename T> friend ostream& operator<<(const MathematicaDebug& S, const T& x)
  {
    if(channel_verbosity[S.channel] < S.level && channel_verbosity["all"] < S.level) return nullstream;

    ostream& stream = channel_stream.find(S.channel) != channel_stream.end()? channel_stream[S.channel]
				      : (channel_stream.find("all") != channel_stream.end()?
					 channel_stream["all"] : cout);
    
    if(prefix>0) stream << "(*" << S.channel << "*) ";
    stream << x;

    return stream;
  }
};

#endif
