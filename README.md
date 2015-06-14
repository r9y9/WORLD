# WORLD

[![Build Status](https://travis-ci.org/r9y9/WORLD.svg?branch=master)](https://travis-ci.org/r9y9/WORLD)
[![Build status](https://ci.appveyor.com/api/projects/status/4j72afijlat4lb8w/branch/master?svg=true)](https://ci.appveyor.com/project/r9y9/world/branch/master)

This repository hosts a slightly modified version of WORLD to provide a easy way to use from external programs. See [here](http://ml.cs.yamanashi.ac.jp/world/english/index.html) for the original WORLD.

This version of WORLD is used by the following bindings:

- [r9y9/WORLD.jl](https://github.com/r9y9/WORLD.jl)
- [r9y9/go-world](https://github.com/r9y9/go-world)
- [sotelo/world.py](https://github.com/sotelo/world.py)

## Supported Platforms

- Linux
- Mac OS X
- Windows

## Installation

```bash
./waf configure && ./waf
sudo ./waf install
```

**NOTE**: MSVC will be used for compilation on windows by default. If you prefer to use `g++`, use `--check-cxx-compiler` option like:

```bash
./waf configure --check-cxx-compiler=g++
```

## Changes from the original WORLD

- Add `extern C` in header files
- Add `DioByOptPtr` (for calling from Julia that doesn't support struct-passing by value)
- Integrate waf
- Support pkg-config
- Add `GetWORLDVersion` function [#4]

[#4]: https://github.com/r9y9/WORLD/pull/4
