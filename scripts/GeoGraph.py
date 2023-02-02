import argparse
import os
import signal
import subprocess
from threading import Timer

BASE_DIR = '/home/funke/devel/geograph'
BUILD_DIR = os.path.join(BASE_DIR, 'build')
EXEC = os.path.join(BUILD_DIR, 'GeoGraph')

DATA_DIR = os.path.join(BASE_DIR, 'data')

DIST_NAME2CHAR = {'uni': 'u', 'gaussian': 'g', 'grid': 'd', 'road': 'r', 'stars': 's', 'circle': 'c', 'bubbles': 'b'}
DIST_CHAR2CHAR = {'u': 'uni', 'g': 'gaussian', 'd': 'grid', 'r': 'road', 's': 'stars', 'c': 'circle', 'b': 'bubbles'}
SEEDS = [8158, 14030, 18545, 20099, 24065, 35700, 37197, 38132, 59135, 60315]

DISTS = ['u', 'g', 'd', 'r', 's', 'c', 'b']
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

    return out, err, proc.returncode
