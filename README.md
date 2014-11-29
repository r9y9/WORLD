# WORLD

This repository hosts a slightly modified version of WORLD to provide a easy way to use from external programs. See [here](http://ml.cs.yamanashi.ac.jp/world/english/index.html) for the original WORLD.

This version of WORLD is used by the following bindings:

- [r9y9/WORLD.jl](https://github.com/r9y9/WORLD.jl)
- [r9y9/go-world](https://github.com/r9y9/go-world)

## Supported Platforms

- Linux
- Mac OS X

Note that the WORLD (probably) works in windows but currently I don't provide any build script. You can build WORLD manually.

## Installation

     ./waf configure && ./waf
     sudo ./waf install

## Changes from the original WORLD

- Add `extern C` in header files
- Change name of old interface `Dio` -> `DioOld`
- Add `DioByOptPtr` (for calling from Julia that doesn't support struct-passing by value)
- Integrate waf
- Support pkg-config
