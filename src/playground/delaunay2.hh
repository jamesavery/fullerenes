class IntrinsicDelaunay: public Triangulation {
public:
  matrix<double> distances;
  matrix<double> edge_lengths_d6y; // zero if !edge
};
