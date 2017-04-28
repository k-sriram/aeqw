#!/usr/bin/env python

from isynspec import INLIN
import subprocess

Z = int(raw_input("Enter the Atomic number: "))
Q = int(raw_input("Enter Ionization: "))
T = float(raw_input("Enter Threshold: "))
MINL = float(raw_input("Enter Min Lambda: "))
MAXL = float(raw_input("Enter Max lambda: "))

lines = []
proc1 = subprocess.Popen('cat synspeclines | grep "\<{0:d}.{1:0>2d}\>"'.format(Z,Q), shell=True, stdout=subprocess.PIPE)
for l in proc1.stdout.read().split('\n'):
    if l:
        lines.append(INLIN(l))
for line in lines:
    if line.GF > T and line.ALAM > MINL and line.ALAM < MAXL:
        print line

#print 'cat synspeclines | grep "\<{0:d}.{1:0>2d}\>"'.format(Z,Q)
