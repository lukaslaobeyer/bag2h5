# bag2h5

Convert rosbag to HDF5.

Written quickly and without much regard for efficiency, but should be
able to process bagfiles on the order of 50 GB within a minute or so.

## Build & run:

Using Nix,
```
nix develop
bazel build //bag2h5
./bazel-bin/bag2h5/bag2h5 in.bag out.h5
```

## Limitations:

* Does not handle variable-length arrays.
* Untested with multiple bag files. Not sure what happens if bags are
  provided out of order.
* Does not use HDF5 compression.
