cc_library(
    name = "highfive",
    hdrs = glob(["include/highfive/**/*.hpp"],),
    strip_include_prefix="include",
    visibility = ["//visibility:public"],
    defines = [
        "H5_USE_EIGEN",
    ],
    deps = [
        "@hdf5",
        "@eigen",
    ],
)
