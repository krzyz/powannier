# POWannier

Simple library for finding exponentially localized expressions of Wannier functions for given potential.

## Prerequisites

* C++ compiler supporting C++14 standard.
* [CMake](https://cmake.org) version 3.7 or later.
* [Armadillo](http://arma.sourceforge.net) linear algebra library.
* [Openmp](https://www.openmp.org) (optional) multithreading support.
* [Doxygen](http://doxygen.nl) (optional) building documentation.
* [Latex](https://www.latex-project.org) (optional) building documentation.

[Cubature](https://github.com/stevengj/cubature) (library for multidimensional integration) and [Catch2](https://github.com/catchorg/Catch2) (testing library) are additional dependencies, these two libraries are however automatically downloaded during CMake build.

## Supported Platforms
 * Linux
 * Microsoft Windows

## Installation

### Resolving dependencies

Armadillo can be installed directly using package managers on some Linux distributions ([Debian](https://packages.debian.org/source/sid/armadillo), [Arch](https://aur.archlinux.org/packages/armadillo/)). In case this option is unavailable the source code can be downloaded from the project site.

When building POWannier using CMake it may be necessary to have built Armadillo library in one of the environment variables checked by CMake. For Armadillo to work, both BLAS and LAPACK libraries must be present in the system (easy to download and install on Linux, there are few versions for Windows, one of them can be found in the Armadillo package).

### Using CMake 

In order to build the library, in the project root directory make and enter build directory:
```
mkdir build && cd build
```

Next configure a project:

```
cmake -DCMAKE_BUILD_TYPE=Release ..
```

Finally, compile the library:

```
make
```

After succesfull execution of `make` command, the built library files can be found in `lib` directory, while `include` directory contains all of the header files.

#### Optional: building the documentation

To build the documentation, `BUILD_DOCUMENTATION` option must be set during configuration:

```
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_DOCUMENTATION=ON ..
```

To finish generating the documentation run:

```
make doc
```

The html files containing the documentation are in `docs/html` directory.

### Running tests

After succesfull compilation, the test binary is located in `bin` directory. To run the tests simply run `./tests` (or `.\tests.exe` on Windows)


## Using the library

Examples can be found in `examples` directory and in the documentation.

## Author

**Krzysztof Biedroń**

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

The method used in this library has been described in the following works:

* Steven Kivelson. *Wannier functions in one-dimensional disordered systems: Application to fractionally charged solitons.* Physical Review B, **26**, 4269, 1982
* Ulf Bissbort. *Dynamical effects and disorder in ultracold bosonic matter.* PhD thesis, 2012.