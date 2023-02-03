#!/usr/bin/python3
import re
import signal
import sys

try:
    import natsort

    fileSorter = natsort.natsorted
except ImportError as e:
    print("python natsort library not installed -- falling back to standard sort")
    fileSorter = sorted

from GeoGraph import *


class tee:
    def __init__(self, _fd1, _fd2):
        self.fd1 = _fd1
        self.fd2 = _fd2

    def __del__(self):
        if self.fd1 != sys.stdout and self.fd1 != sys.stderr:
            self.fd1.close()
        if self.fd2 != sys.stdout and self.fd2 != sys.stderr:
            self.fd2.close()

    def write(self, text):
        self.fd1.write(text)
        self.fd2.write(text)

    def flush(self):
        self.fd1.flush()
        self.fd2.flush()


parser = argparse.ArgumentParser()

# dists to run
parser.add_argument("-d", "--dist",
                    help="point distribution(s) to benchmark[_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar, _c_ircle, _b_ubbles]",
                    choices=DISTS, default='u', nargs='+', type=str)
parser.add_argument("--nMin", help="minimum number of points", type=int)
parser.add_argument("--nMax", help="maximum number of points", type=int)
parser.add_argument("-s", "--seed", help="seed(s) for RNG", nargs='+', type=int)
parser.add_argument("-i", "--its", help="iterations per seed", type=int, default=1)

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
parser.add_argument("-t", "--timelimit", help="time limit for algorithm", type=int, default=1800)  # 30 Minutes
parser.add_argument("", "--stats", help="collect statistics", type=bool, default=False)

args = parser.parse_args()

rePointFilename = re.compile('points_([a-z]+)_([0-9]+)_([0-9]+).csv')

if args.conesMin and args.conesMax:
    args.cones = range(args.conesMin, args.conesMax, 2)
else:
    args.cones = [args.cones]

fTimings = tee(sys.stdout, open("timings.csv", "x"))
print("algorithm kernel k dist n seed rep t failed timeout", file=fTimings)

if args.stats:
    fStats = tee(sys.stdout, open("stats.csv", "x"))
    print("algorithm kernel k dist n seed c steps ipPro isPro delPro maxIpQ maxIsQ maxDelQ maxSlSize failed timeout",
          file=fStats)

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
        lArgs.benchmark = args.its
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

                    timings = []
                    stats = []

                    nPoints = 2147483647
                    timeout = 0
                    failed = 0

                    if rc == 0:
                        # Extract metrics from output
                        for line in out.split('\n'):
                            s = str(line).strip()
                            if "RESULT" in s:
                                nPoints = int(s.split(" n=")[1].split(" ")[0])
                                t = int(s.split(" t=")[1].split(" ")[0])
                                r = int(s.split(" rep=")[1].split(" ")[0])
                                timings.append(t)

                            if "STATS" in s:
                                k = int(s.split(" k=")[1].split(" ")[0])
                                steps = int(s.split(" steps=")[1].split(" ")[0])
                                ipPro = int(s.split(" ipPro=")[1].split(" ")[0])
                                isPro = int(s.split(" isPro=")[1].split(" ")[0])
                                delPro = int(s.split(" delPro=")[1].split(" ")[0])
                                maxIpQ = int(s.split(" maxIpQ=")[1].split(" ")[0])
                                maxIsQ = int(s.split(" maxIsQ=")[1].split(" ")[0])
                                maxDelQ = int(s.split(" maxDelQ=")[1].split(" ")[0])
                                maxSlSize = int(s.split(" maxSlSize=")[1].split(" ")[0])
                                stats.append((k, steps, ipPro, isPro, delPro, maxIpQ, maxIsQ, maxDelQ, maxSlSize))

                    elif rc == -signal.SIGTERM:
                        timings.append(lArgs.timelimit)
                        timeout = 1
                    else:
                        timings.append(lArgs.timelimit)
                        failed = 1

                    # CSV format: algorithm kernel k dist n seed rep t failed timeout
                    for t in enumerate(timings):
                        print(lArgs.alg,
                              lArgs.kernel,
                              lArgs.cones,
                              lArgs.dist,
                              nPoints,
                              lArgs.seed,
                              *t,
                              failed,
                              timeout,
                              sep=" ",
                              file=fTimings)

                    # CSV format: algorithm kernel k dist n seed c steps ipPro isPro delPro maxIpQ maxIsQ maxDelQ maxSlSize failed timeout
                    for s in stats:
                        print(lArgs.alg,
                              lArgs.kernel,
                              lArgs.cones,
                              lArgs.dist,
                              nPoints,
                              lArgs.seed,
                              *s,
                              failed,
                              timeout,
                              sep=" ",
                              file=fStats)
