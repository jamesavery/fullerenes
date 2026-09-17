#pragma once

#include <stdexcept>
#include <string>

// Failure reading or interpreting a geometry file. A catchable std::runtime_error
// subtype thrown by the from_file / from_ply / to_ply / .geo readers and writers so a
// malformed file never takes down the host process. The Code is a closed,
// reason-named set of the modeled failure categories (style-failures.md); .what()
// carries the specifics.
struct mesh_io_error : std::runtime_error {
  enum class Code {
    NullFile,           // no readable / writable stream
    UnsupportedFormat,  // header format unsupported, or a required property/type missing
    MalformedFile,      // truncated body, bad element count, or allocation failure
    InvalidTopology,    // out-of-range index, non-manifold arc, open fan, bad degree, inconsistent winding
    NotATriangulation,  // a valid polyhedron, but the target type requires triangular faces
    FaceTooLarge,       // writer: a face has more than 255 vertices
    UnknownFormat,      // from_file: file extension / format string not recognised
    EmptyGraph,         // from_file: the delegated parser produced zero vertices
    NonSimplicial,      // .geo: a self-loop or parallel edge where a simple graph is required
    IndexOutOfRange,    // .geo: record index >= record count
    CapacityExceeded,   // .geo writer: vertex count or degree outside the file's capacities
    ValueOutOfRange,    // .geo writer: a coordinate that the file's number format cannot hold
    HeaderMismatch,     // .geo append: the options differ from the file's header
  };
  Code code;
  mesh_io_error(Code code_, const std::string &what_) : std::runtime_error(what_), code(code_) {}
};
