#!/usr/bin/env python
# autoeqw.py
# Program to automate the calculation of abundance of an element given its
# equivalent width in a particular spectral line.
# K.Sriram
# Created: 20/04/2017

from isynspec import *
import math




INFN = 'aeqw.in'
OUTFN = 'aeqw.out'
INITABUN = 1e-4
NULLABUN = 1e-10
RANGE = 5.0
EPSILON = 0.1
IS = ISynspec()


allLines = [] # allLines stores the information of all the lines that are going to be used
# testLines stores which all lines have to be tested. It is a list. Each
# element can either be a string or a or a testline. A testLine has the 
# format: Couple of a tuple and a float. The tuple is a list of lines to be 
# tested simultaneously. The float is the target equivalent width.
testLines = []
tempTL = []

lineNo = 1

with open(INFN) as f:
    line = f.readline()
    tokens = line.split()
# The first line has the format:
# modelFN INLTE SPACE VTB BROAD LOGHE
    IS.modelFN = tokens[0]
    IS.INLTE = int(tokens[1])
    IS.SPACE = float(tokens[2])
    IS.VTB = float(tokens[3])
    BROAD = float(tokens[4])
    LOGHE = float(tokens[5])
    for line in f:
        lineNo += 1
        if line[0] == '#': # Unprinted Comment
            continue
        if line[0] == 'C': # Printed Comment
            testLines.append(line[2:])
            continue
        tokens = line.split()
        if not tokens: # Empty Line
            continue
        if len(tokens) < 11 or len(tokens) > 12:
            raise InvalidInput('\nLine {0:d}: {1}'.format(lineNo,line))
        try:
            if len(tokens) == 11:
                allLines.append(INLIN(line))
            else:
                allLines.append(INLIN(' '.join(tokens[:-1])))
                tempTL.append(allLines[-1])
                if float(tokens[-1]) != 0:
                    testLines.append( (tempTL,float(tokens[-1])) )
                    tempTL = []
        except InvalidInput as e:
            print e
            raise InvalidInput('\nLine {0:d}: {1}'.format(lineNo,line))
        except Exception:
            print "Error while processing line {0:d}\n{1}".format(lineNo,line)
            raise

INCONSISTENT = False
if not IS.TestTG():
    print ''
    print("         *************************** ERROR ***************************")
    print("         *Discrepancy between fort.8 and model temperature / gravity.*")
    print("         *************************** ERROR ***************************")
    print ''
    INCONSISTENT = True
allLines.sort()
IS.LINELIST = allLines
IS.write19()
# Determining the bounds of the synthetic spectrum ALAM0 and ALAM1. The multiplication by ten is for conversion from nm to A. Also writing to 19 and 55

def InitParam(testLine, allLines):
    global IS
    IS.ALAM0 = min(testLine).ALAM * 10 - RANGE
    IS.ALAM1 = max(testLine).ALAM * 10 + RANGE
    
    IS.write55()
    
def Overlap(bin,box):
    l = max(bin[0],box[0])
    r = min(bin[1],box[1])
    if l > r:
        return 0.
    return (r-l)/(bin[1] - bin[0])


# Calculate the Equivalent width of a particular line
def CalcEqw(testLine):
    global IS
    box = ( min(testLine).ALAM * 10 - BROAD, max(testLine).ALAM * 10 + BROAD )
    total = 0
    for bin in IS.EQW:
        total += bin[1] * Overlap(bin[0],box)
    return total

# Set the abundance and run SYNSPEC and read the output
def Run(abundances):
    global IS
    IS.ABUNDANCES = abundances
    IS.write56()
    IS.run()
    IS.read16()
finAbun = [] 
# Iterating over all testLines
for tl in testLines:
    if type(tl) == str:
        finAbun.append('')
        continue
    testLine = tl[0]
    xeqw = tl[1]
    Z = testLine[0].Z
    print 'Calculating for following lines with target equivalent width:', xeqw
    for t in testLine:
        print str(t)[:-1]
    InitParam(testLine, allLines)
    # Setting Zero
    Run([(Z,NULLABUN)])
    zero = CalcEqw(testLine)

    # Finding the abundance that gives reasonable eqw
    trials = [INITABUN]
    results = []
    while not results or abs(results[-1] - xeqw) > EPSILON:
        Run([(Z,trials[-1])])
        results.append(CalcEqw(testLine) - zero)
        if results[-1] == 0:
            if trials[-1] < .1:
                trials.append(trials[-1] * 10)
                continue
            else:
                finAbun.append('No line Detected')
                break
        elif results[-1] * xeqw < 0:
            if trials[-1] > 0.1:
                finAbun.append('Line Strength Insufficient')
            else:
                finAbun.append('Emmision/Absorption mismatch')
            break
        else:
            trials.append(xeqw * trials[-1] / results[-1])
    else:
        finAbun.append('{0:0.2e} {1:.2f}'.format(trials[-1], math.log(trials[-1],10) + LOGHE))
    print 'Result:', finAbun[-1], '\n'

# Writing the output
with open(OUTFN,'w') as f:
    if INCONSISTENT:
        f.write("Model Inconsistent ")
    f.write("{0:.2f} {1:.2f}\n".format(IS.TEMP,IS.LOGG))
    f.write("LAMBDANM  Z.Q  ABUN/He  LOGABUN")
    for i in range(len(testLines)):
        if type(testLines[i]) == str:
            f.write('\n' + testLines[i].rstrip('\n'))
            continue
        for line in testLines[i][0]:
            f.write('\n{0:.4f} {1:2d}.{2:0>2d}'.format(line.ALAM,line.Z,line.Q))
        f.write(' {0}'.format(finAbun[i]))
    f.write('\n')
print "Total runs:", IS.runs
