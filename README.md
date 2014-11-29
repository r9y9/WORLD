# WORLD
------------------------------------

This repository hosts a slightly modified version of WORLD to provide a easy way to use from external programs. See [here](http://ml.cs.yamanashi.ac.jp/world/english/index.html) for the original WORLD.

## Changes from the original WORLD

- Add `extern C` in header files
- Change name of old interface `Dio` -> `DioOld`
- Add `DioByOptPtr` (for calling from Julia that doesn't support struct-passing by value)
- Integrate waf
- Support pkg-config

## Supported Platforms

- Linux
- Mac OS X

Note that the WORLD (probably) works in windows but currently I don't provide any build script. You can build WORLD manually.

## Installation

     ./waf configure && ./waf
     sudo ./waf install

## Bindings

- [Julia](https://github.com/r9y9/WORLD.jl)
- [Golang](https://github.com/r9y9/go-world)