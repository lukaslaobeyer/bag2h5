cc_binary(
    name = "bag2h5",
    srcs = [
    	"type_info.h",
        "bag2h5.cc"
    ],
    deps = [
        ":logging",
        "@eigen",
        "@magic_enum",
        "@embag",
        "@highfive",
    ],
)

cc_library(
    name = "logging",
    srcs = ["logging.cc"],
    hdrs = ["logging.h"],
)
