#!/usr/bin/env python3
"""LO cross section (pb) of a point from the quirkgen logs, averaged over seeds. usage: xsec.py <dir> <tag-glob>"""
import glob, re, sys
v = [float(re.search(r"sigma = (\S+) pb", open(f).read()).group(1)) for f in glob.glob("%s/%s.gen.log" % (sys.argv[1], sys.argv[2]))]
print(sum(v) / len(v))
