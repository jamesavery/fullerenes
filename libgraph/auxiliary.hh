#ifndef AUXILIARY_HH
# define AUXILIARY_HH
#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <sstream>
#include <assert.h>
using namespace std;

#define insert_unique(v,x) if(std::find(v.begin(),v.end(),x) == v.end()) v.push_back(x); 

template <typename S, typename T> ostream& operator<<(ostream& s, const pair<S,T>& p)
{
  s << "{" << p.first << "," << p.second << "}";
  return s;
}

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

template<typename K, typename V> vector<V> get_values(const vector<pair<K,V> >& m)
{
  vector<V> values(m.size());
  int i=0;
  for(typename vector<pair<K,V> >::const_iterator kv(m.begin()); kv!=m.end(); kv++,i++)
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

template <typename T> T from_string(const string& s)
{
  stringstream S(s);
  T x;
  S >> x;
  return x;
}

string pad_string(const string& s, int length, char padchar = '0');

string filename_extension(const string& filename);


int gcd(int a, int b);

template <typename T> vector<T> operator*(const vector<T>& xs, const T& x)
{
  vector<T> ys(xs.size());
  for(int i=0;i<xs.size();i++) ys[i] = xs[i] * x;
  return ys;
}

template <typename T> vector<T> operator+(const vector<T>& xs, const T& x)
{
  vector<T> ys(xs.size());
  for(int i=0;i<xs.size();i++) ys[i] = xs[i] + x;
  return ys;
}

template <typename T> vector< vector<T> > operator+(const vector< vector<T> >& xs, const T& x)
{
  vector< vector<T> > ys(xs.size());
  for(int i=0;i<xs.size();i++) ys[i] = xs[i] + x;
  return ys;
}


template <typename T> class IDCounter: public map<T,int> {
public:
  int nextid;
  vector<T> reverse;

  IDCounter(int start=0) : nextid(start) {}
  
  int insert(const T& x){
    typename map<T,int>::const_iterator it(map<T,int>::find(x));
    if(it != this->end()) return it->second;
    else {
      map<T,int>::insert(make_pair(x,nextid));
      reverse.push_back(x);
      return nextid++;
    }
  }
  
  const T& invert(int idx) const { 
    assert(idx>=0 && idx<nextid);
    return reverse[idx];
  }

  int operator()(const T& x) const {
    typename map<T,int>::const_iterator it(map<T,int>::find(x));
    if(it != this->end()) return it->second;
    else return -1;
  }
};


#endif
