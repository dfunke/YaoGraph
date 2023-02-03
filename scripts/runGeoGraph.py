#!/usr/bin/python3

from GeoGraph import *

parser = argparse.ArgumentParser()

# generate points
parser.add_argument("-d", "--dist",
                    help="point distribution [_u_ni, _g_aussian, gri_d_, _r_oad, _s_tar, _c_ircle, _b_ubbles]",
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
