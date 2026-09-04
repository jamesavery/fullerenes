#pragma once
// Fortran formatted I/O edit descriptors, emulated column-exactly.
//
// Legacy data files (the isomer database, the program's listings) are
// written and read by Fortran FORMAT statements; this is the vocabulary
// for reading and writing such records from C++, one edit descriptor at a
// time:  Aw (characters), Iw (integer), Fw.d (fixed-point real).
//
// Input follows Fortran BLANK='NULL' semantics: blanks anywhere in a
// numeric field are ignored and an all-blank field reads as 0; a sign must
// precede the digits.  Fw.d input accepts only a plain decimal with one
// point (Fortran would also take an exponent, or apply the implied scale
// to a field without a point -- neither occurs in files written by an F
// edit descriptor, and the implied scale would be a silent misread).
// Output right-justifies the value in the field; a value that does not fit
// throws (Fortran prints asterisks and silently loses it), and Fw.d output
// reproduces gfortran: the optional leading zero of |v| < 1 is dropped
// when that is what it takes to fit ("-.50000" in F7.5).
//
// Each reader takes the record s and the current column pos (0-based),
// consumes len columns and advances pos.
// @throws std::invalid_argument on a malformed or missing input field, or
//         output that does not fit.
#include <string>
#include <string_view>

namespace fortran {

// The next len columns of s; the descriptor names the edit descriptor in
// the diagnostic.
std::string_view field(const std::string& s, int& pos, int len, const char* descriptor);

void        read_A(char* result, const std::string& s, int& pos, int len);
long        read_I(const std::string& s, int& pos, int len);
double      read_F(const std::string& s, int& pos, int len);
std::string write_I(long v, int len);
std::string write_F(double v, int len, int decimals);

}  // namespace fortran
