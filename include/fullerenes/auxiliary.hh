#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <sstream>
#include <algorithm>
#include <sys/stat.h>
#include <assert.h>
#include <unordered_map>
#include <unordered_set>

extern char LIST_OPEN;
extern char LIST_CLOSE;

using namespace std;

typedef int node_t;
typedef vector< vector<node_t> > neighbours_t;
typedef vector< bool > edges_t;

template <typename S, typename T> ostream& operator<<(ostream& s, const pair<S,T>& p);

#define container_output(container) \
  template <typename T> ostream& operator<<(ostream& s, const container<T>& v) \
  { \
  s << LIST_OPEN; \
  for(typename container<T>::const_iterator x(v.begin());x!=v.end();){ \
    s << *x; \
    if(++x!=v.end()) s << ","; \
  } \
  s << LIST_CLOSE; \
  return s; \
}

container_output(vector);
container_output(list);
container_output(set);


// Directed edge is an ordered pair of nodes
typedef pair<node_t,node_t> arc_t;

struct edge_t : public pair<node_t,node_t> {
  edge_t() {}
  edge_t(const pair<node_t,node_t>& p) : pair<node_t,node_t>(std::min(p.first,p.second),std::max(p.first,p.second)) {}
  edge_t(const node_t u, const node_t v): pair<node_t,node_t>(std::min(u,v),std::max(u,v)) {}
  edge_t(const int index) {
    node_t u=0;
    for(;u*(u-1)/2<=index;u++) ;
    u--;
    first = u;
    second = index-u*(u-1)/2;
  }
  inline size_t index() const { 
    const node_t v = first, u = second;
    return u*(u-1)/2 + v; 
  }
};



#define insert_unique(v,x) if(std::find(v.begin(),v.end(),x) == v.end()) v.push_back(x); 

template <typename S, typename T> ostream& operator<<(ostream& s, const pair<S,T>& p)
{
  s << LIST_OPEN << p.first << "," << p.second << LIST_CLOSE;
  return s;
}


// TODO: Macro instead of repeating
template<typename K, typename V> vector<K> get_keys(const map<K,V>& m)
{
  vector<K> keys(m.size());
  int i=0;
  for(const auto &kv: m)
    keys[i++] = kv.first;
  return keys;
}

template<typename K, typename V> vector<K> get_keys(const unordered_map<K,V>& m)
{
  vector<K> keys(m.size());
  int i=0;
  for(const auto &kv: m)
    keys[i++] = kv.first;
  return keys;
}

template<typename K, typename V> vector<K> get_keys(const vector<pair<K,V> >& m)
{
  vector<K> keys(m.size());
  int i=0;
  for(const auto &kv: m)
    keys[i++] = kv.first;
  return keys;
}


template<typename K, typename V> vector<V> get_values(const map<K,V>& m)
{
  vector<V> values(m.size());
  int i=0;
  for(const auto &kv: m)
    values[i++] = kv.second;
  return values;
}

template<typename K, typename V> vector<V> get_values(const unordered_map<K,V>& m)
{
  vector<V> values(m.size());
  int i=0;
  for(const auto &kv: m)
    values[i++] = kv.second;
  return values;
}

template<typename K, typename V> vector<V> get_values(const vector<pair<K,V> >& m)
{
  vector<V> values(m.size());
  int i=0;
  for(const auto &kv: m)
    values[i++] = kv.second;
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
  for(size_t i=0;i<xs.size();i++) ys[i] = xs[i] + x;
  return ys;
}

template <typename T> vector<T> operator-(const vector<T>& xs, const T& x)
{
  vector<T> ys(xs.size());
  for(int i=0;i<xs.size();i++) ys[i] = xs[i] - x;
  return ys;
}

template <typename T> vector<T> operator+(const vector<T>& xs, const vector<T>& ys)
{
  vector<T> result(xs.size());
  for(int i=0;i<xs.size();i++) result[i] = xs[i] + ys[i];
  return result;
}

template <typename T> vector<T> operator-(const vector<T>& xs, const vector<T>& ys)
{
  vector<T> result(xs.size());
  for(int i=0;i<xs.size();i++) result[i] = xs[i] - ys[i];
  return result;
}

template <typename T> vector<T> operator-(const vector<T>& xs)
{
  vector<T> result(xs.size());
  for(int i=0;i<xs.size();i++) result[i] = -xs[i];
  return result;
}

template <typename T> vector<T> sorted(const vector<T>& xs)
{
  vector<T> result(xs.begin(), xs.end());
  std::sort(result.begin(), result.end()); 
  return result;
}

template <typename T> vector< vector<T> > operator+(const vector< vector<T> >& xs, const T& x)
{
  vector< vector<T> > ys(xs.size());
  for(size_t i=0;i<xs.size();i++) ys[i] = xs[i] + x;
  return ys;
}

template <typename T> T sum(const vector<T>& xs)
{
  T sum = 0;
  for(const auto &x: xs) sum += x;
  return sum;
}

template <typename T> T mean(const vector<T>& xs)
{
  T sum_ = sum(xs);
  return sum_/xs.size();
}

template<typename T> void hash_combine(size_t &seed, T const &key) {
  hash<T> hasher;
  seed ^= hasher(key) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std {
  template<typename T1, typename T2> struct hash<pair<T1, T2>> {
    size_t operator()(const pair<T1, T2> &p) const {
      size_t seed(0);
      hash_combine(seed, p.first);
      hash_combine(seed, p.second);
      return seed;
    }
  };

  template<typename IntType> struct hash<vector<IntType>> { // Vectors of integers smaller than 32 bit
    size_t operator()(const vector<IntType> &v) const {
      return std::hash<u32string>()(u32string(v.begin(),v.end()));      
    }
  };

}

template <typename T> class IDCounter: public unordered_map<T,size_t> {
public:
  size_t nextid;
  vector<T> reverse;

  IDCounter(size_t start=0) : nextid(start) {}
  IDCounter(const vector<T>& xs, size_t start=0) : nextid(start+xs.size()), reverse(xs)
  {
    for(size_t i=0;i<xs.size();i++) unordered_map<T,size_t>::insert({xs[i],i+start});
  }
  
  size_t insert(const T& x){
    typename unordered_map<T,size_t>::const_iterator it(unordered_map<T,size_t>::find(x));
    if(it != this->end()) return it->second;
    else {
      unordered_map<T,size_t>::insert({x,nextid});
      reverse.push_back(x);
      return nextid++;
    }
  }
  
  const T& invert(size_t idx) const { 
    assert(idx>=0 && idx<nextid);
    return reverse[idx];
  }

  size_t operator()(const T& x) const {
    typename unordered_map<T,size_t>::const_iterator it(unordered_map<T,size_t>::find(x));
    if(it != this->end()) return it->second;
    else return -1;
  }
};


// C++-style getline with posix files.
bool getline(FILE *file, string& str);



string trim(const string& str, const string& wschars);

template <typename T> vector<T> split(const string& parse_str, const string& delimiters, const string wschars=" \t\r\n")
{
  vector<string> string_result = split<string>(parse_str,delimiters,wschars);
  vector<T> result;

  for(string s: string_result) result.push_back(from_string<T>(s));

  return result;
}

// Version of split that handles empty strings ("a;;b;c;" -> {"a","","b","c",""} in stead of {"a","b","c"}.
template <> vector<string> split(const string& parse_str, const string& delimiters, const string wschars);

// Fast double-ended queue implemented with cyclic buffer and possibility of look-ahead and look-back
template<typename T> class Deque: public std::vector<T> {
public:
  const size_t N;
  size_t front_index, back_index;
  bool   empty() const { return front_index==back_index; }
  size_t size()  const { return back_index-front_index;  }
  void   clear()       { front_index=back_index; }
  
  const T& front(size_t offset=0) const {
    assert(!empty() && (offset<size()));
    return (*this)[front_index+offset];
  }
  const T& back (size_t offset=0) const {
    assert(!empty() && (offset<size()));
    return (*this)[(back_index-1-offset+N)%N];
  }
  
  T pop_front() {
    assert(!empty());
    T x((*this)[front_index]);
    front_index = (front_index+1) % N;
    return x;
  }
  T pop_back()  {
    assert(!empty());
    T x((*this)[(back_index-1+N)%N]);
    back_index  = (back_index-1+N)% N;
    return x;
  }

  void push_front(const T& x) {
    (*this)[front_index] = x; front_index = (front_index-1+N) % N;
    assert(!empty());		// Overflow-check
  }
  void push_back (const T& x) {
    (*this)[back_index]  = x; back_index  = (back_index+1)    % N;
    assert(!empty());		// Overflow-check
  }

  Deque(size_t max_length) : std::vector<T>(max_length), N(max_length), front_index(0), back_index(0) {}

  friend ostream& operator<<(ostream &s, const Deque q)
  {
    vector<T> contents(&q[q.front_index],&q[q.back_index]);
    return (s << contents);
  }
};

size_t file_size(const char *filename);
size_t file_size(const string filename);
