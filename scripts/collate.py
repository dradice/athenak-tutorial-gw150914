#!/usr/bin/env python

import argparse
from numpy import inf
import sys

parser = argparse.ArgumentParser()
parser.add_argument("input", metavar="input.dat", nargs='+',
        help="Input files to collate")
parser.add_argument("-c", "--comment", dest="comment", action="append",
        default=['#'], metavar='C', help="Comment character (default: \'#\')")
parser.add_argument("-e", "--epsilon", dest="epsilon", type=float,
        default=1e-15, help="Precision used in comparing timestamps "
        "(default: 1e-15)")
parser.add_argument("--header", dest="header", default=0,
        type=int, metavar='N', help='Number of header lines to include')
parser.add_argument("--tmax", dest="tmax", default=inf, type=float,
        help="Maximum time to include (default: inf)")
parser.add_argument("--include-all-comments", dest="include_all_comments",
        action="store_true", help="Include all comments in the output file")
parser.add_argument("-i", "--include-comments", dest="include_comments",
        action='store_true', help="Include comments from the first file only")
parser.add_argument("-o", "--output", dest="output", required=True,
        help="Output file", metavar="output.dat")
parser.add_argument("-t", "--time-idx", dest="tidx", type=int,
        default=1, help="Index of the time column, from 1 (default: 1)")
args = parser.parse_args()
args.comment.insert(0, '\n')

# Insert header
ifile = open(args.input[0], 'r')
ofile = open(args.output, 'w')
for i in range(args.header):
    ofile.write(ifile.readline())
ifile.close()

# Collate data
told  = None
for fnum, fname in enumerate(args.input):
    for dline in open(fname, 'r'):
        skip = False
        for c in args.comment:
            if dline[:len(c)] == c:
                if args.include_all_comments or (
                        args.include_comments and fnum == 0):
                    ofile.write(dline)
                skip = True
                break
        if skip:
            continue

        tnew = float(dline.split()[args.tidx-1])
        if tnew > args.tmax:
            break
        if told is None or tnew > told*(1 + args.epsilon):
            ofile.write(dline)
            told = tnew
