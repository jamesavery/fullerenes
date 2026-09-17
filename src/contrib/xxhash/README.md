# xxHash 0.8.1 (vendored)

`xxhash.h` is the single-header xxHash library by Yann Collet, version 0.8.1
(`XXH_VERSION_MAJOR/MINOR/RELEASE` = 0/8/1), BSD 2-Clause licensed — see the
notice at the top of the file. Source: https://github.com/Cyan4973/xxHash;
this copy was taken unmodified from the vulkan-validationlayers-dev package
(md5 63854394ee5e31db5681712c6b2527d0).

Used by `src/c++/geo-format.cc` for the `.geo` checksum (GEO-FORMAT.md §9),
compiled with `XXH_INLINE_ALL`, so no xxHash symbol is exported. The XXH3
output is frozen since 0.8.0; the checksum depends on that.
