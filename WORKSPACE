load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")
load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")

# zlib
http_archive(
    name = "zlib",
    build_file = "@//deps/zlib:BUILD",
    sha256 = "b3a24de97a8fdbc835b9833169501030b8977031bcb54b3b3ac13740f846ab30",
    strip_prefix = "zlib-1.2.13",
    urls = [
        "https://zlib.net/zlib-1.2.13.tar.gz",
        "https://mirror.bazel.build/zlib.net/zlib-1.2.13.tar.gz",
    ],
)

# Eigen
http_archive(
    name = "eigen",
    build_file = "@//deps/eigen:BUILD",
    sha256 = "8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72",
    url = "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz",
    strip_prefix="eigen-3.4.0"
)

# Magic Enum
http_archive(
    name = "magic_enum",
    sha256 = "62bd7034bbbfc3d7806001767d5775ab42f3ff33bb38366e1ceb21102f0dff9a",
    url = "https://github.com/Neargye/magic_enum/archive/refs/tags/v0.8.2.tar.gz",
    strip_prefix="magic_enum-0.8.2"
)

# HDF5
http_archive(
    name = "hdf5",
    url = "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.0/src/hdf5-1.14.0.tar.gz",
    sha256 = "a571cc83efda62e1a51a0a912dd916d01895801c5025af91669484a1575a6ef4",
    strip_prefix = "hdf5-1.14.0",
    build_file = "@//deps/hdf5:BUILD",
    patches = [ "@//deps/hdf5:disable_subfiling.patch" ],
    patch_args = [ "-p1" ],
)
http_archive(
    name = "highfive",
    sha256 = "ab51b9fbb49e877dd1aa7b53b4b26875f41e4e0b8ee0fc2f1d735e0d1e43d708",
    build_file = "@//deps/highfive:BUILD",
    url = "https://github.com/BlueBrain/HighFive/archive/v2.6.2.tar.gz",
    strip_prefix="HighFive-2.6.2",
)

# Embag and dependencies:
git_repository(
    name = "com_github_nelhage_rules_boost",
    commit = "4dcbe9e03c20da67195aef61a76f73a0b30b9bf9",
    remote = "https://github.com/nelhage/rules_boost",
)
load("@com_github_nelhage_rules_boost//:boost/boost.bzl", "boost_deps")
boost_deps()

http_archive(
    name = "embag",
    sha256 = "394ed9cb5cee14b7472e2f763403bcecb59197bbe3ad492058e9a815eec43c06",
    build_file = "@//deps/embag:BUILD",
    url = "https://github.com/embarktrucks/embag/archive/74c0b5f9d50bd45bcb6ed8e44718cd60924c13d0.tar.gz",
    strip_prefix="embag-74c0b5f9d50bd45bcb6ed8e44718cd60924c13d0/lib",
    patches = [
        "@//deps/embag:add_raw_data_accessor.patch",
        "@//deps/embag:use_std_span.patch"
    ],
    patch_args = [ "-p1" ],
)

http_archive(
    name = "liblz4",
    build_file = "@//deps/lz4:BUILD",
    sha256 = "658ba6191fa44c92280d4aa2c271b0f4fbc0e34d249578dd05e50e76d0e5efcc",
    strip_prefix = "lz4-1.9.2",
    urls = ["https://github.com/lz4/lz4/archive/v1.9.2.tar.gz"],
)

http_archive(
    name = "libbz2",
    build_file = "@//deps/bz2:BUILD",
    sha256 = "ab5a03176ee106d3f0fa90e381da478ddae405918153cca248e682cd0c4a2269",
    strip_prefix = "bzip2-1.0.8",
    urls = ["https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz"],
)
