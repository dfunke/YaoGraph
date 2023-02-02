#!/usr/bin/python3
import re
import signal

try:
    import natsort
    fileSorter = natsort.natsorted
except ImportError as e:
    print("python natsort library not installed -- falling back to standard sort")
    fileSorter = sorted

from GeoGraph import *

parser = argparse.ArgumentParser()

# dists to run
parser.add_argument("-d", "--dist",
                    help="point distribution(s) to benchmark[_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar, _c_ircle, _b_ubbles]",
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

for dDir in fileSorter(os.listdir(DATA_DIR)):

    if not dDir in DIST_NAME2CHAR or not DIST_NAME2CHAR[dDir] in args.dist:
        continue

    for pFile in fileSorter(os.listdir(os.path.join(DATA_DIR, dDir))):
        match = rePointFilename.match(pFile)

        lArgs = argparse.Namespace

        lArgs.infile = os.path.join(DATA_DIR, dDir, pFile)

        lArgs.dist = DIST_NAME2CHAR[match.group(1)]
        lArgs.n = int(match.group(2))
        lArgs.seed = int(match.group(3))
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

                    out, err, rc = run_algorithm(lArgs)

                    print(out)

                    if err:
                        print(err)

                    if rc == -signal.SIGTERM:
                        print("Time Out")
                    else:
                        print("Failed")
