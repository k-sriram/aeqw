#!/usr/bin/env python
# autoeqw.py
# Program to automate the calculation of abundance of an element given its
# equivalent width in a particular spectral line.
# K.Sriram
# Created: 20/04/2017

from time import time
import logging
import logging.handlers
from isynspec import *
import math

startTime = time()


# Initializing the logger
logger = logging.getLogger('aeqw')
logger.setLevel(logging.DEBUG)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
console.setFormatter(logging.Formatter('%(levelname)-8s: %(message)s'))
logger.addHandler(console)
filelog = logging.handlers.RotatingFileHandler('aeqw.log','a',maxBytes=0,backupCount=9)
filelog.setLevel(logging.DEBUG)
filelog.setFormatter(logging.Formatter('%(asctime)s - %(name)-15s %(levelname)-8s: %(message)s'))
logger.addHandler(filelog)
filelog.doRollover()

logger.debug('Running program: Automatic Equation width solver.')
logger.debug('Initialization')

INFN = 'aeqw.in'
OUTFN = 'aeqw.out'
INITABUN = 1e-4
NULLABUN = 1e-10

logger.debug(' Parameters: INFN: {0:s}, OUTFN, {1:s}, INITABUN: {2:.1e}, NULLABUN: {3:.1e}'.format(INFN,OUTFN,INITABUN,NULLABUN))

BROAD = 2.0
RANGE = 5.0
EPSILON = 0.1
IS = ISynspec()

INITABUNZWISE = {}
IS.read56()
for  i in IS.ABUNDANCES:
    INITABUNZWISE[i[0]] = i[1]

allLines = [] # allLines stores the information of all the lines that are going to be used
# testLines stores which all lines have to be tested. It is a list. Each
# element can either be a string or a or a testline. A testLine has the 
# format: Couple of a tuple and a float. The tuple is a list of lines to be 
# tested simultaneously. The float is the target equivalent width.
testLines = []
tempTL = []

lineNo = 1

logger.debug("Reading input file.")
with open(INFN) as f:
    line = f.readline()
    tokens = line.split()
# The first line has the format:
# modelFN INLTE SPACE VTB BROAD LOGHE
    IS.modelFN = tokens[0]
    LOGHE = float(tokens[1])
    BROAD = float(tokens[2])
    RANGE = float(tokens[3])
    EPSILON = float(tokens[4])
    logger.debug(' Parameters from input file: Model filename: {0}'.format(IS.modelFN))
    logger.debug(' > LOG[He]: {0:f}, BROAD: {1:f}, RANGE: {2:f}, EPSILON: {3:e}'.format(LOGHE,BROAD,RANGE,EPSILON))
    for line in f:
        logger.debug(' Processing line: {0}'.format(line.strip()))
        lineNo += 1
        if line[0] == '#': # Unprinted Comment
            continue
        if line[0] == 'C': # Printed Comment
            testLines.append(line[2:])
            continue
        if line.strip() == '': # Empty Line
            continue
        try:
            inline = INLIN(line)
            allLines.append(inline)
            if inline.remainder.strip():
                tempTL.append(allLines[-1])
                if float(inline.remainder) != 0:
                    testLines.append( (tempTL,float(inline.remainder)) )
                    tempTL = []
        except Exception:
            logger.error("Error while processing line {0:d}\n{1}\n".format(lineNo,line))
            raise

INCONSISTENT = False
if not IS.TestTG():
    logger.warning("Discrepancy between fort.8 and model temperature / gravity.")
    INCONSISTENT = True
allLines.sort()
IS.LINELIST = allLines
IS.write19()
# Determining the bounds of the synthetic spectrum ALAM0 and ALAM1. The multiplication by ten is for conversion from nm to A. Also writing to 19 and 55

def InitParam(testLine, allLines):
    global IS
    IS.ALAM0 = min(testLine).ALAM * 10 - RANGE
    IS.ALAM1 = max(testLine).ALAM * 10 + RANGE
    logger.debug(" InitParam: Setting range of synthetic spectrum: ({0:.1f}, {1:.1f})".format(IS.ALAM0,IS.ALAM1))
    
    IS.write55()
    
def Overlap(bin,box):
    l = max(bin[0],box[0])
    r = min(bin[1],box[1])
    if l > r:
        return 0.
    return math.sqrt((r-l)/(bin[1] - bin[0]))
def Secant(x,f,y):
    if abs(f[-1] - f[-2]) < EPSILON:
        return -1
    return (x[-2]*(f[-1]-y) - x[-1]*(f[-2]-y))/(f[-1]-f[-2])

# Calculate the Equivalent width of a particular line
def CalcEqw(testLine):
    global IS
    if len(IS.EQW) < 2:
        logger.warning("  CalcEqw: SYNSPEC did not generate output in fort.16")
        return None
    box = ( min(testLine).ALAM * 10 - BROAD, max(testLine).ALAM * 10 + BROAD )
    logger.debug("  CalcEqw: Calculating Equivalent width; including bins in {0}.".format(str(box)))
    total = 0
    for bin in IS.EQW:
        total += bin[1] * Overlap(bin[0],box)
    logger.debug("  CalcEqw: eqw = {0:f}".format(total))
    return total

# Set the abundance and run SYNSPEC and read the output
def Run(abundances):
    global IS
    logger.debug("  Setting abundance: {0}".format(str(abundances)))
    IS.ABUNDANCES = abundances
    IS.write56()
    IS.run()
    IS.read16()
finAbun = [] 
# Iterating over all testLines
logger.debug("Estimating abundance for all lines.")

for tl in testLines:
    if type(tl) == str:
        finAbun.append('')
        continue
    testLine = tl[0]
    xeqw = tl[1]
    Z = testLine[0].Z
    logger.info('Calculating for following lines with target equivalent width: %f', xeqw)
    for t in testLine:
        logger.info(str(t))
    InitParam(testLine, allLines)
    # Setting Zero
    logger.debug(" Performing zero check")
    abun = NULLABUN
    Run([(Z,abun)])
    zero = CalcEqw(testLine)
    while zero is None:     # If the program didn't compute the bins
        IS.RELOP /= 10
        logger.debug(" > Setting RELOP parameter to {0:.1e}".format(IS.RELOP))
        IS.write55()
        Run([(Z,abun)])
        zero = CalcEqw(testLine)
    logger.debug(" > Zero = {0:f}".format(zero))

    # Finding the abundance that gives reasonable eqw
    trials = [INITABUNZWISE.get(Z,INITABUN)]
    results = []
    while not results or abs(results[-1] - xeqw) > EPSILON:
        logger.debug(" Running for abundance: {0:e}, target width: {1:f}".format(trials[-1],xeqw))
        Run([(Z,trials[-1])])
        results.append(CalcEqw(testLine) - zero)
        if results[-1] is None or results[-1] == 0:
            if trials[-1] < .1:
                logger.debug("Negligible width detected, multiplying abundance by 10")
                trials.append(trials[-1] * 10)
                continue
            else:
                finAbun.append('No line Detected')
                logger.debug('No line Detected')
                break
        elif results[-1] * xeqw < 0:
            logger.debug('eqw * xeqw < 0')
            if trials[-1] > 0.1:
                finAbun.append('Line Strength Insufficient')
            else:
                finAbun.append('Emmision/Absorption mismatch')
            break
        else:
            logger.debug("  Guess = {0:e}, Result = {1:f}, Target = {2:f}, Diff = {3:f}, Epsilon = {4:f}".format(trials[-1],results[-1],xeqw,xeqw-results[-1],EPSILON))
            if len(results) < 2:
                trials.append(xeqw * trials[-1] / results[-1])
                logger.debug(" Using linear approximation for new guess: {0:e}".format(trials[-1]))
            else:
                trials.append(Secant(trials,results,xeqw))
                if trials[-1] < 0:
                    trials[-1] = (xeqw * trials[-2] / results[-1])
                    logger.debug(" Using linear approximation for new guess: {0:e}".format(trials[-1]))
                else:
                    logger.debug(" Using secant method for new guess: {0:e}".format(trials[-1]))
    else:
        finAbun.append('{0:0.2e} {1:.2f}'.format(trials[-1], math.log(trials[-1],10) + LOGHE))
    logger.info('Result: %s', finAbun[-1])

# Writing the output
with open(OUTFN,'w') as f:
    logger.debug("Writing Output")
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
logger.info("Total runs: %d", IS.runs)

# Resetting the fort.56 file
IS.ABUNDANCES = []
for i in INITABUNZWISE:
    IS.ABUNDANCES.append((i,INITABUNZWISE[i]))
IS.write56()

logger.info("Runtime: {0:.3f}".format(time()-startTime))
