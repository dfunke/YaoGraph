#!/usr/bin/python3
import argparse
import os
import os.path
import signal
import subprocess
import re
from threading import Timer

BASE_DIR = '/home/funke/devel/geograph'
BUILD_DIR = os.path.join(BASE_DIR, 'build')
EXEC = os.path.join(BUILD_DIR, 'GeoGraph')

DATA_DIR = os.path.join(BASE_DIR, 'data')

DIST_NAME2CHAR = {'uni': 'u', 'gaussian': 'g', 'grid': 'd', 'road': 'r', 'stars': 's'}
DIST_CHAR2CHAR = {'u': 'uni', 'g': 'gaussian', 'd': 'grid', 'r': 'road', 's': 'stars'}
SEEDS = [8158, 14030, 18545, 20099, 24065, 35700, 37197, 38132, 59135, 60315]

DISTS = ['u', 'g', 'd', 'r', 's']
ALGS = ['s', 'g', 'n', 'c']
KERNS = ['i', 'c', 'p']


def run_algorithm(args: argparse.Namespace):
    params = [EXEC,
              "--alg", str(args.alg),
              "--kernel", args.kernel,
              "--cones", str(args.cones),
              "--cellOcc", str(args.cellOcc),
              ]

    if args.infile:
        params.extend(["--infile", args.infile])
        if args.dist:
            params.extend(["--dist", args.dist])  # for naming
        if args.seed:
            params.extend(["--seed", str(args.seed)])  # for naming

    elif args.n:
        params.extend(["--dist", args.dist])
        params.extend(["--n", str(args.n)])
        params.extend(["--seed", str(args.seed)])

    if args.outfile:
        params.extend(["--outfile", args.outfile])

    if args.stdout:
        params.extend(["--stdout"])

    if args.benchmark:
        params.extend(["--benchmark"])

    proc = subprocess.Popen(params, stdout=subprocess.PIPE, universal_newlines=True, preexec_fn=os.setsid)

    def kill_proc():
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)

    t = Timer(args.timelimit, kill_proc)
    t.start()
    out, err = proc.communicate()
    t.cancel()

    return out, err


parser = argparse.ArgumentParser()
# benchmark
parser.add_argument("-b", "--benchmark", help="run benchmark suite", action='store_true')

# generate points
parser.add_argument("-d", "--dist", help="point distribution [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar]",
                    choices=DISTS, default='u', type=str)
parser.add_argument("-n", "--n", help="number of points to generate", type=int)
parser.add_argument("-s", "--seed", help="seed for RNG", default=8158, type=int)

# points file
parser.add_argument("-f", "--infile", help="file with points", type=str);

# algorithm to use
parser.add_argument("-a", "--alg", help="algorithm to use [_s_weepline, _g_rid, _n_aive, _c_gal]",
                    choices=ALGS, default='s', type=str)
parser.add_argument("-c", "--kernel",
                    help="kernel to use [_i_nexact, CGALExact_p_redicatesInexactConstructions, CGALExactpredicatesInexact_c_onstructions]",
                    choices=KERNS, default='i', type=str)
parser.add_argument("-k", "--cones", help="number of cones, default 6", default=6, type=int)
parser.add_argument("-g", "--cellOcc", help="number of points per cell (grid algorithm)", default=100, type=int)

# output
parser.add_argument("-o", "--outfile", help="file for graph output", type=str)
parser.add_argument("-p", "--stdout", help="write graph to stdout", action='store_true')

# execution
parser.add_argument("-t", "--timelimit", type=int, default=1800)  # 30 Minutes

args = parser.parse_args()

out, err = run_algorithm(args)

print(out)

if err:
    print(err)
