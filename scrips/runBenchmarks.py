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

# dists to run
parser.add_argument("-d", "--dist",
                    help="point distribution(s) to benchmark[_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar]",
                    choices=DISTS, default='u', nargs='+', type=str)
parser.add_argument("--nMin", help="minimum number of points", type=int)
parser.add_argument("--nMax", help="maximum number of points", type=int)
parser.add_argument("-s", "--seed", help="seed(s) for RNG", nargs='+', type=int)

# algorithm to use
parser.add_argument("-a", "--alg", help="algorithm(s) to use [_s_weepline, _g_rid, _n_aive, _c_gal]",
                    choices=ALGS, default='s', nargs='+', type=str)
parser.add_argument("-c", "--kernel",
                    help="kernel(s) to use [_i_nexact, CGALExact_p_redicatesInexactConstructions, CGALExactpredicatesInexact_c_onstructions]",
                    choices=KERNS, default='i', nargs='+', type=str)
parser.add_argument("-k", "--cones", help="number of cones, default 6", default=6, type=int)
parser.add_argument("--conesMin", help="minimum number of cones", type=int)
parser.add_argument("--conesMax", help="maximum number of cones", type=int)
parser.add_argument("-g", "--cellOcc", help="number of points per cell (grid algorithm)", default=100, type=int)

# execution
parser.add_argument("-t", "--timelimit", type=int, default=1800)  # 30 Minutes

args = parser.parse_args()

rePointFilename = re.compile('points_([a-z]+)_([0-9]+)_([0-9]+).csv')

if args.conesMin and args.conesMax:
    args.cones = range(args.conesMin, args.conesMax, 2)
else:
    args.cones = [args.cones]

for dDir in os.listdir(DATA_DIR):

    if not DIST_NAME2CHAR[dDir] in args.dist:
        continue

    for pFile in os.listdir(os.path.join(DATA_DIR, dDir)):
        match = rePointFilename.match(pFile)

        lArgs = argparse.Namespace

        lArgs.infile = pFile

        lArgs.dist = DIST_NAME2CHAR[match.group(1)]
        lArgs.n = match.group(2)
        lArgs.seed = match.group(3)
        lArgs.cellOcc = args.cellOcc
        lArgs.outfile = None
        lArgs.stdout = None
        lArgs.benchmark = True
        lArgs.timelimit = args.timelimit

        if args.nMin and lArgs.n < args.nMin:
            continue
        if args.nMax and lArgs.n > args.nMax:
            continue

        if args.seed and not lArgs.seed in args.seed:
            continue

        for alg in args.alg:
            for kern in args.kernel:
                for con in args.cones:

                    lArgs.alg = alg
                    lArgs.kernel = kern
                    lArgs.cones = con

                    out, err = run_algorithm(lArgs)
                    print(out)

                    if err:
                        print(err)
