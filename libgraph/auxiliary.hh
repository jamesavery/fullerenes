#ifndef AUXILIARY_HH
# define AUXILIARY_HH
#include <iostream>
using namespace std;

// TODO: There are a number of functions in this file that are not particularly
// pertaining to geometry, but are just "miscellaneous" stuff. Move to a more
// fitting place.
#define container_output(container) \
  template <typename T> ostream& operator<<(ostream& s, const container<T>& v) \
  { \
  s << "{"; \
  for(typename container<T>::const_iterator x(v.begin());x!=v.end();){ \
    s << *x; \
    if(++x!=v.end()) s << ","; \
  } \
  s << "}"; \
  return s; \
}




container_output(vector);
container_output(list);
container_output(set);


template <typename S, typename T> ostream& operator<<(ostream& s, const pair<S,T>& p)
{
  s << "{" << p.first << "," << p.second << "}";
  return s;
}

template<typename K, typename V> vector<K> get_keys(const map<K,V>& m)
{
  vector<K> keys(m.size());
  int i=0;
  for(typename map<K,V>::const_iterator kv(m.begin()); kv!=m.end(); kv++,i++)
    keys[i] = kv->first;
  return keys;
}

template<typename K, typename V> vector<K> get_keys(const vector<pair<K,V> >& m)
{
  vector<K> keys(m.size());
  int i=0;
  for(typename vector<pair<K,V> >::const_iterator kv(m.begin()); kv!=m.end(); kv++,i++)
    keys[i] = kv->first;
  return keys;
}


template<typename K, typename V> vector<V> get_values(const map<K,V>& m)
{
  vector<V> values(m.size());
  int i=0;
  for(typename map<K,V>::const_iterator kv(m.begin()); kv!=m.end(); kv++,i++)
    values[i] = kv->second;
  return values;
}

template<typename from, typename to> vector<to> convert_vector(const vector<from>& x)
{
  vector<to> y(x.size());
  for(int i=0;i<x.size();i++) y[i] = x[i];
  return y;
}

template<typename K, typename V> K getFirst(const pair<K,V>& x){ return x.first; }
template<typename K, typename V> V getSecond(const pair<K,V>& x){ return x.second; }

template<typename S, typename T> pair<T,S> reverse(const pair<S,T>& p){ return pair<T,S>(p.second,p.first); }

template <typename T> int sgn(const T& val) { return (T(0) < val) - (val < T(0)); }

// Undirected edge is an unordered pair of nodes

template <typename T> string to_string(const T& x)
{
  ostringstream s;
  s << x;
  return s.str();
}

string pad_string(const string& s, int length, char padchar = '0');

#endif
