# Yao Graph Generator

This tool generates the Yao-Graph of a given point set.
The Yao graph connects each point to its nearest neighbor in each of _k_ cones drawn around it.
This implementation works for two-dimensional inputs.

## Compilation

YaoGraph requires a modern C++ compiler with C++20 support.
To use exact predicates and constructions the CGAL library is required.

```shell
mkdir build && cd build
cmake ../
make
```

To use CGAL's kernels for predicates and constructions turn on the `WITH_CGAL` CMake option

```shell
mkdir build && cd build
cmake -DWITH_CMAKE=ON ../
make
```

Further CMake options are:
- `WITH_CAIRO` to draw images of the generated Yao graphs
- `WITH_STATS` to collect statistics during the sweepline execution (performance penalty)
- `WITH_TESTS` compile gtest-Suite
- `ROAD_DIR` directory containing DIMACS road networks (see below)
- `STAR_DIR` directory containing Gaia star catalogue (see below)

## Usage

### YaoGraph

YaoGraph can either randomly generate points or use an input file with points

Usage
```shell
./YaoGraph --help

Yao Graph generator:
  -h, --help                produce help message
  -b, --benchmark arg (=1)  run algorithm i times
  -d, --dist arg (=u)       point distribution [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar, _c_ircle, _b_ubbles]
  -n, --n arg               number of points to generate
  -s, --seed arg (=8158)    seed for RNG
  -f, --infile arg          file with points
  -a, --alg arg (=s)        algorithm to use [_s_weepline, _g_rid, _n_aive, _c_gal]
  -c, --kernel arg (=i)     kernel to use [_i_nexact, CGALExact_p_redicatesInexactConstructions, CGALExactpredicatesInexact_c_onstructions]
  -k, --cones arg (=6)      number of cones, default 6
  -g, --cellOcc arg (=100)  number of points per cell (grid algorithm)
  -o, --outfile arg         file for graph output
  -p, --stdout              write graph to stdout
  -v, --verify              verify graph (-v grid-based, -v -v naive yao (slow))
  --stats                   collect execution statistics (requires WITH_STATS compile flag)

```

Input files have the following format
```shell
# n 1000
# b 0 0 1 1
0.227937 0.682662
0.396135 0.946997
0.844108 0.681706
0.728227 0.537687
```

The first line specifies the number of points in the file, 
the second line defines the bounding box of the contained points.

### Generator

The generator application can be used to randomly generate input points and store the resulting points in files.
Especially the road and star distribution require some time to generate (see Data)

```shell
./Generator --help

Point Generator:
  -h, --help              produce help message
  -a, --all               generate all available distributions
  -n, --n arg             number of points to generate
  --minN arg (=1000)      minimum number of points to generate
  --maxN arg (=10000000)  maxium number of points to generate
  -s, --seed arg (=8158)  seed for RNG
  -i, --inst arg          number of instances

list of distributions to generate [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar, _c_ircle, _b_ubbles]
```

## Benchmarks

### Data

All synthetic distributions can be generated without further inputs.
- For road networks you need to download the desired network from http://www.diag.uniroma1.it//challenge9/download.shtml
- The Gaia 2 star data can be downloaded from https://www.cosmos.esa.int/web/gaia/data-release-2

Benchmarks can then be executed with the Python script `scripts/runBenchmarks.py`:

```shell
usage: runBenchmarks.py [-h] [-d {u,g,d,r,s,c,b} [{u,g,d,r,s,c,b} ...]] [--nMin NMIN] [--nMax NMAX]
                        [-s SEED [SEED ...]] [-i ITS] [-a {s,g,n,c} [{s,g,n,c} ...]] [-c {i,c,p} [{i,c,p} ...]]
                        [-k CONES] [--conesMin CONESMIN] [--conesMax CONESMAX] [-g CELLOCC] [-t TIMELIMIT]
                        [--stats]

options:
  -h, --help            show this help message and exit
  -d {u,g,d,r,s,c,b} [{u,g,d,r,s,c,b} ...], --dist {u,g,d,r,s,c,b} [{u,g,d,r,s,c,b} ...]
                        point distribution(s) to benchmark[_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar, _c_ircle,
                        _b_ubbles]
  --nMin NMIN           minimum number of points
  --nMax NMAX           maximum number of points
  -s SEED [SEED ...], --seed SEED [SEED ...]
                        seed(s) for RNG
  -i ITS, --its ITS     iterations per seed
  -a {s,g,n,c} [{s,g,n,c} ...], --alg {s,g,n,c} [{s,g,n,c} ...]
                        algorithm(s) to use [_s_weepline, _g_rid, _n_aive, _c_gal]
  -c {i,c,p} [{i,c,p} ...], --kernel {i,c,p} [{i,c,p} ...]
                        kernel(s) to use [_i_nexact, CGALExact_p_redicatesInexactConstructions,
                        CGALExactpredicatesInexact_c_onstructions]
  -k CONES, --cones CONES
                        number of cones, default 6
  --conesMin CONESMIN   minimum number of cones
  --conesMax CONESMAX   maximum number of cones
  -g CELLOCC, --cellOcc CELLOCC
                        number of points per cell (grid algorithm)
  -t TIMELIMIT, --timelimit TIMELIMIT
                        time limit for algorithm
  --stats               collect statistics

```

### Plotting

Plots can be generated with `benchmark/plot.py`.




