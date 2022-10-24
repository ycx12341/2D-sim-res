# 2D-Sim (C++ Implementation)

// TODO

## Usage

### Executable Sample

// TODO

### Compile the code

- Required C++ Standard: **17**
  - [Compiler support for C++17](https://en.cppreference.com/w/cpp/compiler_support/17)
- May consider to provide a `Makefile`.

> <details>
> <summary><b>Ubuntu: Update your GCC/G++ compiler</b></summary>
> 
> 1. Adding PPA source:
>    `add-apt-repository ppa:ubuntu-toolchain-r/test`
> 2. Update apt:
>    `sudo apt update`
> 3. Install GCC/G++ v11:
>    `sudo apt install gcc-11 g++-11`
> 4. Switch to installed GCC/G++ version:
>    `sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 10`
>    `sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 10`
> 5. Check version number:
>    `gcc --version`
>    `g++ --version`
> 
> </details>

#### Unix

```
<g++> main.cpp sample/scc.cpp seed/nrom.cpp seed/sample.cpp seed/seed.cpp -o <executable_file_name> -lpthread
```

#### Windows

```
<g++> main.cpp sample/scc.cpp seed/nrom.cpp seed/sample.cpp seed/seed.cpp -o <executable_file_name>
```