// geo_io.hh -- opening .geo files by path for the bindings (GEO-FORMAT.md). The C++
// API takes FILE*; these open and close the stream and report an unopenable path
// as OSError (FileNotFoundError, PermissionError, ...).
#pragma once

#include <pybind11/pybind11.h>

#include <cerrno>
#include <cstdio>
#include <memory>
#include <string>

#include <fcntl.h>
#include <unistd.h>

namespace pyf {

namespace py = pybind11;

struct file_closer { void operator()(FILE* f) const { std::fclose(f); } };
using unique_file = std::unique_ptr<FILE, file_closer>;

[[noreturn]] inline void raise_os_error(const std::string& path, int err) {
    errno = err;
    PyErr_SetFromErrnoWithFilename(PyExc_OSError, path.c_str());
    throw py::error_already_set();
}

inline unique_file open_for_reading(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) raise_os_error(path, errno);
    return unique_file(f);
}

// Run `write` on a stream for `path` and close it; a failed close is a failed write.
// Neither mode truncates at open: a fresh write (geo::write) truncates only once
// every record has been encoded, so a refused write leaves an existing file as it
// was -- and removes a file this call created. Appending opens read/write and
// creates the file when absent, never truncating another writer's file.
template<class F>
bool write_file(const std::string& path, bool append, F&& write) {
    int fd = ::open(path.c_str(), (append ? O_RDWR : O_WRONLY) | O_CREAT | O_EXCL, 0666);
    const bool created = fd >= 0;
    if (!created && errno == EEXIST) fd = ::open(path.c_str(), append ? O_RDWR : O_WRONLY);
    if (fd < 0) raise_os_error(path, errno);
    FILE* f = ::fdopen(fd, append ? "r+b" : "wb");
    if (!f) { const int err = errno; ::close(fd); raise_os_error(path, err); }
    unique_file file(f);
    bool ok;
    try {
        ok = write(file.get());
    } catch (...) {
        file.reset();
        if (created && !append) ::unlink(path.c_str());
        throw;
    }
    return std::fclose(file.release()) == 0 && ok;
}

}  // namespace pyf
