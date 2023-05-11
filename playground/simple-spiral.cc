#include "fullerenes/fullerenegraph.hh"
#include "fullerenes/triangulation.hh"
#include "fullerenes/polyhedron.hh"
#include "fullerenes/spiral.hh"

#define INLINE __attribute__((always_inline)) inline

// Statically allocated deque with sub-byte sized data type.
// NB: Caller is in charge of checking that a full deque is never pushed,
//     and an empty queue is never popped. Either by calling full() resp. empty()
//     before pushing resp. popping, or by ensuring structurally that it Can Never Happen(tm).
template <int nbits> struct packed_deque {
private:
  static constexpr uint8_t mask[9] = {0x0,0x1,0x3,0x7,0xf,0x1f,0x3f,0x7f,0xff};
  ssize_t start, end, capacity, length;
  uint8_t *data;

public:
  // data must be preallocated memory of (capacity*nbits)/8+1 bytes.
  packed_deque(int capacity, uint8_t *data) : start(0), end(0),
					      capacity(capacity), length(0), data(data) {
    for(ssize_t i=0;i<capacity*nbits/8+1;i++) data[i] = 0;
  }

  INLINE ssize_t size() const { return length; }
  INLINE bool  full()   const { return length==capacity; }
  INLINE bool  empty()  const { return length==0; }

  INLINE static ssize_t bits_index(ssize_t i){ return (i*nbits)>>3;  }
  INLINE static ssize_t bits_start(ssize_t i){ return (i*nbits)&0x7; }

  // read_element(i) reads the i'th n-bit field in the flat memory, i.e. without reference to start and end.
  INLINE uint8_t read_element(const ssize_t i) const {
    constexpr uint8_t m      = mask[nbits];
    const size_t ix          = bits_index(i), shift = bits_start(i); // Byte index in memory, and bit index for start of value
    const uint16_t *word_ptr = (uint16_t*)(data+ix);          // We need to load two bytes, as value may cross a byte boundary

    uint16_t word = *word_ptr;
    return (word >> shift) & m;
  }

  // write_element(i,x) write x to the i'th n-bit field in the flat memory, i.e. without reference to start and end.  
  INLINE void write_element(const ssize_t i, const uint16_t x) {
    constexpr uint16_t m = mask[nbits];

    const size_t ix    = bits_index(i), shift = bits_start(i); // Byte index in memory, and bit index for start of value

    uint16_t *word_ptr = (uint16_t*)(data+ix);          // We need to load two bytes, as value may cross a byte boundary
    uint16_t word = *word_ptr;   

    word &= ~(m<<shift);     // Zero out old field
    word |= (x&m)<<shift;    // Write x into field
    *word_ptr = word;	     // Store.
  }

  INLINE void push_back(const uint8_t& x) {
    if(length>0) end = (end+1)%capacity; // Whenever length<=1, we want front() == back()
    
    write_element(end,x);
    length++;
  }

  INLINE void push_front(const uint8_t& x) {
    if(length>0) start = (start+capacity-1)%capacity; // Whenever length<=1, we want front() == back()

    write_element(start,x);
    length++;
  }

  INLINE uint8_t front() const { return read_element(start); }  
  INLINE uint8_t back()  const { return read_element(end);   }

  INLINE uint8_t pop_front() {
    uint8_t x = front();
    length--;
    start = (start+1)%capacity;	
    return x;
  }

  INLINE uint8_t pop_back() {
    uint8_t x = back();
    length--;
    end = (end+capacity-1) % capacity;
    return x;
  }

  friend ostream &operator<<(ostream& S, const packed_deque& D)
  {
    S << "{capacity:"<<D.capacity<<", length:"<<D.length
      <<", start:"<<D.start<<", end:"<<D.end<<"}";
    S<< "[";
    ssize_t data_length = ((D.length*nbits)/8)+(((D.length*nbits)%8)!=0);

    for(ssize_t i=0;i<D.length;i++)
      S << int(D.read_element((D.start+i)%D.capacity)) << (i+1<D.length?",":"");
    S << "]";
    return S;
  }
};


bool spiral(const Triangulation& g, tri_t start_tri)
{
  uint8_t boundary_code_data[(g.N*3)/8+1];
  packed_deque<3> boundary_code(g.N,boundary_code_data);

  cout << "0: boundary_code = " << boundary_code << ";\n";
  boundary_code.push_back(3);
  cout << "1: boundary_code = " << boundary_code << ";\n";  
  boundary_code.push_back(4);
  cout << "2: boundary_code = " << boundary_code << ";\n";
  boundary_code.push_front(2);
  cout << "3: boundary_code = " << boundary_code << ";\n";    
  boundary_code.push_front(1);
  cout << "3: boundary_code = " << boundary_code << ";\n";  
  boundary_code.push_front(0);
  cout << "4: boundary_code = " << boundary_code << ";\n";

  //  while(!boundary_code.empty()){
  for(int i=0;i<200;i++){
    uint8_t x;
    x = boundary_code.pop_back();
    boundary_code.push_front(x);
    cout << "boundary_code = " << boundary_code << ";\n";
  }
  cout << "boundary_code = " << boundary_code << ";\n";
  
  return true;
}

int main(int ac, char **av)
{
  if(ac<2){
    fprintf(stderr,"Syntax: %s \"<spiral-name>\n"
	    "Example: %s C100-[1,2,3,4,5,6,47,48,49,50,51,52]-fullerene\n\n",av[0],av[0]);
    return -1;
  }
  
  spiral_nomenclature name(av[1]);

  Triangulation g(name);

  node_t u = 0, v = g.neighbours[u][0], w = g.next_on_face(u,v);
  spiral(g,{u,v,w});
  //  cout << "g = " << g << ";\n";
}
