# The post-day-3 simulation in the second experiment with SCC invasion patterns (C++ Implementation)

Another option to boost the computational efficiency is to reimplement the algorithm into a language faster than R. Taking the post-day 3 simulation in our second experiment with SCC invasion patterns, the computational efficiency was boosted by around $84.2\%$ in our C++ version. 

Welcome to extending the code to simulate invasion patterns on other days (or other datasets).

## Usage

### Executable Sample

See [Release](https://github.com/ycx12341/2D-sim-res/releases).

### Compile the code

- Required C++ Standard: **17**
  - [Compiler support for C++17](https://en.cppreference.com/w/cpp/compiler_support/17)
- May consider to provide a `Makefile`.

> <details>
> <summary><b>[FCQ] Ubuntu: Update your GCC/G++ compiler</b></summary>
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

## Code Structure

```
\algo				Main algorithms
	calculate.h		ESS; BW; ABC_BCD; SSE
	pars.h			Structures for parameters to be evaluated
	pattern.h		Simulation pattern
	pde.h			PDE machanism
	scc.h			Simulation for SCC
	sim2d.h
\collection
	collection.h	Extended functions for collections
	matrix.h		Highly customized matrix designed for Sim2D
\red_den			Cell density
\sample				Sampling 
	sample.h		
	scc.cpp			
\seed				Randomness
	nrom.cpp		Normal Distribution
	sample.cpp		Random Samples, Permutations and Uniform Distribution
	seed.cpp		Seeding and Randomness.
	seed.h			
config.h
```

