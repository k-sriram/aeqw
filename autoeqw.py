#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# autoeqw.py
# Program to automate the calculation of abundance of an element given its
# equivalent width in a particular spectral line.
# K.Sriram
# Created: 20/04/2017

version = "1.2.dev1"

import sys
from time import time
from datetime import datetime
import logging
import logging.handlers
import math
from configparser import ConfigParser, ExtendedInterpolation
from argparse import ArgumentParser
import json
from isynspec import ISynspec, INLIN, aeqwISError

CONFFN = "aeqw.conf"
confgetter = {"str": "get", "float": "getfloat", "int": "getint", "bool": "getboolean"}


conf = ConfigParser(interpolation=ExtendedInterpolation())
conf["DEFAULT"] = {
    "INFN": "aeqw.in",
    "OUTFN": "aeqw.out",
    "EXTRALOGFN": "",
    "INITABUN": 1e-4,
    "NULLABUN": 1e-10,
    "LOGATREF": 11.54,
    "BROAD": 2.0,
    "RANGE": 5.0,
    "EPSILON": 0.1,
    "SEP19": False,
    "OUTFMT": "txt",
}
conf["TYPES"] = {
    "INFN": "str",
    "OUTFN": "str",
    "EXTRALOGFN": "str",
    "INITABUN": "float",
    "NULLABUN": "float",
    "LOGATREF": "float",
    "BROAD": "float",
    "RANGE": "float",
    "EPSILON": "float",
    "SEP19": "bool",
    "OUTFMT": "str",
}
conf["aeqw"] = {}

argparser = ArgumentParser(
    description="Program to automate equivalent width finding using SYNSPEC."
)
argparser.add_argument(
    "model",
    help="Name of the input model (such as 'hhe35lt').",
    nargs="?",
    default="fort",
)
argparser.add_argument("-i", "--infn", help="Custom input filename.")
argparser.add_argument("-o", "--outfn", help="Custom output filename.")
argparser.add_argument(
    "-c", action="append", help="Extra configuration options files. Can be repeated."
)
argparser.add_argument("-l", "--extralogfn", help="Extra log file.")
argparser.add_argument(
    "--initabun",
    type=float,
    help="Assumed initial abundance for elements if not described in 'fort.56'.",
)
argparser.add_argument(
    "--nullabun",
    type=float,
    help="Abundance of elements used to estimate equivalent width at zero abundance.",
)
argparser.add_argument("--logatref", type=float, help="Abundance of reference element.")
argparser.add_argument(
    "--broad",
    type=float,
    help="Half width (in Å) upto which absorption is assumed to come from the line. Set it so as to cover the entire line. If set correctly 'wing%%' should be low for isolated lines.",
)
argparser.add_argument(
    "--range",
    type=float,
    help="Half width (in Å) of the generated synthetic spectrum used for analysis.",
)
argparser.add_argument(
    "--epsilon",
    type=float,
    help="Accuracy to which the program will try to match the equivalent width.",
)
argparser.add_argument(
    "--sep19",
    action="store_true",
    help="Whether to generate a separate 'fort.19' file for every line group.",
)
argparser.add_argument(
    "--outfmt",
    help="Format of the output file, Valid options are txt and json.",
    choices=("txt", "json"),
)
argparser.add_argument("-v", action="store_true", help="Print version number and exit.")
args = argparser.parse_args()

argoptions = (
    "infn",
    "outfn",
    "extralogfn",
    "initabun",
    "nullabun",
    "broad",
    "range",
    "epsilon",
    "outfmt",
)

startTime = time()


def readconf(conffn, conf=conf, warnifnotfound=False):
    if conf.read(conffn) == []:
        if warnifnotfound:
            logger.warn(" Configuration file {0:s} not found.".format(conffn))
        else:
            logger.info(
                " Configuration information can be added to file {0:s}".format(conffn)
            )
        return
    logger.info(" Reading Configuration file: {0:s}".format(conffn))


def getconf(param, conf=conf, sec="aeqw"):
    if param in conf[sec]:
        try:
            return getattr(conf[sec], confgetter[conf["TYPES"][param]])(param)
        except ValueError:
            logger.error("Invalid type for parameter {} in configuration".format(param))
            raise
    else:
        logger.error("Parameter {} not found".format(param))
        raise ValueError


def Overlap(bin, box):
    l = max(bin[0], box[0])
    r = min(bin[1], box[1])
    if l > r:
        return 0.0
    return math.sqrt((r - l) / (bin[1] - bin[0]))


def Secant(x, f, y):
    if abs(f[-1] - f[-2]) < getconf("EPSILON"):
        return -1
    return (x[-2] * (f[-1] - y) - x[-1] * (f[-2] - y)) / (f[-1] - f[-2])


# Output formatters


def outputtxt(outputData, outfn):
    with open(outfn, "w") as f:
        f.write("{model:s} {temperature:.2f} {logg:.2f}\n".format_map(outputData[0]))
        if "unit55" in outputData[0]:
            f.write(
                "".join(
                    [
                        "{} = {}\n".format(key, value)
                        for key, value in outputData[0]["unit55"]
                    ]
                )
            )
        f.write("LAMBDANM   Z.Q      Teqw  ABUN/ref  LOGABUN   wing%")
        for row in outputData[1]:
            if row["type"] == "comment":
                f.write("\n" + row["value"])
            elif row["type"] == "line":
                f.write(
                    "".join(
                        [
                            "\n{wavelength: >8.4f}  {ion}".format_map(line)
                            for line in row["lines"]
                        ]
                    )
                )
                abuntxt = "{relabun: >8.2e}  {logabun: >7.2f}   {wingpercent: >4.0f}%".format_map(
                    row["abundance"]
                )
                f.write(" {0:8.2f}  {1}".format(row["target"], abuntxt))
        f.write("\n")


def outputjson(outputData, outfn):
    txt = json.dumps(outputData)
    with open(outfn, "w") as f:
        f.write(txt)


outputformatter = {
    "txt": outputtxt,
    "json": outputjson,
}

# Initializing the logger
logger = logging.getLogger("aeqw")
logger.setLevel(logging.DEBUG)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
console.setFormatter(logging.Formatter("%(levelname)-8s: %(message)s"))
logger.addHandler(console)
filelog = logging.handlers.RotatingFileHandler(
    "aeqw.log", "a", maxBytes=0, backupCount=9
)
filelog.setLevel(logging.DEBUG)
filelog.setFormatter(
    logging.Formatter("%(asctime)s - %(name)-15s %(levelname)-8s: %(message)s")
)
logger.addHandler(filelog)
filelog.doRollover()

logger.info("Running aeqw version: {}".format(version))

if args.v == True:
    logger.debug("Exiting because only version information was asked for.")
    sys.exit(0)

logger.info(str(datetime.now()))

readconf(CONFFN, conf)

if args.c is not None:
    for c in args.c:
        readconf(c, conf, True)

extralog = (
    args.extralogfn
    if args.extralogfn is not None
    else (getconf("EXTRALOGFN") if getconf("EXTRALOGFN") != "" else None)
)

if extralog is not None:
    xtrafilelog = logging.FileHandler(extralog)
    xtrafilelog.setLevel(logging.DEBUG)
    xtrafilelog.setFormatter(
        logging.Formatter("%(asctime)s - %(name)-15s %(levelname)-8s: %(message)s")
    )
    logger.addHandler(xtrafilelog)

for option in argoptions:
    if getattr(args, option) is not None:
        conf["aeqw"][option.upper()] = getattr(args, option)
if args.sep19 == True:
    conf["aeqw"]["SEP19"] = True

logger.info("Running program: Automatic Equation width solver.")
logger.info("Model: {}".format(args.model))
logger.debug("Initialization")

logger.debug(" Parameters:")
for param in conf["aeqw"].keys():
    logger.debug("  {0:s} : {1:s}".format(param, str(getconf(param))))

try:
    with ISynspec(args.model) as IS:

        if "unit55" in conf:
            for param in conf["unit55"]:
                if hasattr(IS, param.upper()):
                    setattr(
                        IS,
                        param.upper(),
                        type(getattr(IS, param.upper()))(conf["unit55"][param]),
                    )
                    logger.info(
                        "Setting unit 55 parameter {} to {}".format(
                            param.upper(), getattr(IS, param.upper())
                        )
                    )

        allLines = (
            []
        )  # allLines stores the information of all the lines that are going to be used
        # testLines stores which all lines have to be tested. It is a list. Each
        # element can either be a string or a or a testline. A testLine has the
        # format: Couple of a tuple and a float. The tuple is a list of lines to be
        # tested simultaneously. The float is the target equivalent width.
        testLines = []
        tempTL = []

        lineNo = 1

        logger.debug("Reading input file.")
        with open(getconf("INFN")) as f:
            for line in f:
                logger.debug(" Processing line: {0}".format(line.strip()))
                lineNo += 1
                if line[0] == "#":  # Unprinted Comment
                    continue
                if line[0] == "C":  # Printed Comment
                    testLines.append(line[2:])
                    continue
                if line.strip() == "":  # Empty Line
                    continue
                try:
                    inline = INLIN(line)
                    allLines.append(inline)
                    if inline.remainder.strip():
                        tempTL.append(allLines[-1])
                        if float(inline.remainder) != 0:
                            testLines.append((tempTL, float(inline.remainder)))
                            tempTL = []
                except Exception:
                    logger.error(
                        "Error while processing line {0:d}\n{1}\n".format(lineNo, line)
                    )
                    raise

        INCONSISTENT = False
        allLines.sort(key=lambda x: x.ALAM)
        IS.LINELIST = allLines
        IS.write19()
        # Determining the bounds of the synthetic spectrum ALAM0 and ALAM1. The multiplication by ten is for conversion from nm to A. Also writing to 19 and 55

        def InitParam(testLine, allLines):
            IS.ALAM0 = min([line.ALAM for line in testLine]) * 10 - getconf("RANGE")
            IS.ALAM1 = max([line.ALAM for line in testLine]) * 10 + getconf("RANGE")
            logger.debug(
                " InitParam: Setting range of synthetic spectrum: ({0:.1f}, {1:.1f})".format(
                    IS.ALAM0, IS.ALAM1
                )
            )

            IS.write55()

            if getconf("SEP19") == True:
                IS.LINELIST = testLine
                IS.write19()

        # Calculate the Equivalent width of a particular line
        def CalcEqw(testLine):
            if len(IS.EQW) < 2:
                logger.warning("  CalcEqw: SYNSPEC did not generate output in fort.16")
                return None, 0
            box = (
                min([line.ALAM for line in testLine]) * 10 - getconf("BROAD"),
                max([line.ALAM for line in testLine]) * 10 + getconf("BROAD"),
            )
            logger.debug(
                "  CalcEqw: Calculating Equivalent width; including bins in {0}.".format(
                    str(box)
                )
            )
            total = 0
            alltotal = 0
            for bin in IS.EQW:
                total += bin[1] * Overlap(bin[0], box)
                alltotal += bin[1]
            logger.debug(
                "  CalcEqw: eqw = {0:f}, alleqw = {1:f}".format(total, alltotal)
            )
            return total, alltotal

        # Set the abundance and run SYNSPEC and read the output
        def Run(abundances):
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
                finAbun.append({"result": "comment"})
                continue
            testLine = tl[0]
            xeqw = tl[1]
            Z = testLine[0].Z
            logger.info(
                "Calculating for following lines with target equivalent width: %f", xeqw
            )
            for t in testLine:
                logger.info(str(t))
            InitParam(testLine, allLines)
            # Setting Zero
            logger.debug(" Performing zero check")
            abun = getconf("NULLABUN")
            Run([(Z, abun)])
            zero, allzero = CalcEqw(testLine)
            while zero is None:  # If the program didn't compute the bins
                if IS.RELOP > 1e-12:
                    IS.RELOP /= 10
                    logger.debug(
                        " > Setting RELOP parameter to {0:.1e}".format(IS.RELOP)
                    )
                    IS.write55()
                Run([(Z, abun)])
                zero, allzero = CalcEqw(testLine)

            logger.debug(" > Zero = {0:f}, allZero = {1:f}".format(zero, allzero))

            # Finding the abundance that gives reasonable eqw
            trials = [IS.INITABUNZWISE.get(Z, getconf("INITABUN"))]
            results = []
            while not results or abs(results[-1] - xeqw) > getconf("EPSILON"):
                logger.debug(
                    " Running for abundance: {0:e}, target width: {1:f}".format(
                        trials[-1], xeqw
                    )
                )
                Run([(Z, trials[-1])])
                results.append(CalcEqw(testLine)[0] - zero)
                if results[-1] is None or (
                    results[-1] >= 0 and results[-1] < xeqw / 10
                ):  # Limiting our increase by a factor of ten
                    if trials[-1] < 0.1:
                        logger.debug(
                            "Negligible width detected, multiplying abundance by 10"
                        )
                        trials.append(trials[-1] * 10)
                        continue
                    else:
                        finAbun.append(
                            {"result": "error", "message": "No line Detected."}
                        )
                        logger.warning("No line Detected")
                        break
                elif results[-1] * xeqw < 0:
                    logger.warning("eqw * xeqw < 0")
                    if trials[-1] > 0.1:
                        finAbun.append(
                            {
                                "result": "error",
                                "message": "Line Strength Insufficient with Emmision/Absorbtion mismatch",
                            }
                        )
                    else:
                        finAbun.append(
                            {
                                "result": "error",
                                "message": "Emmision/Absorption mismatch",
                            }
                        )
                    break
                else:
                    logger.debug(
                        "  Guess = {0:e}, Result = {1:f}, Target = {2:f}, Diff = {3:f}, Epsilon = {4:f}".format(
                            trials[-1],
                            results[-1],
                            xeqw,
                            xeqw - results[-1],
                            getconf("EPSILON"),
                        )
                    )
                    if len(results) < 2 or all(
                        [
                            not (i is None or (i >= 0 and i < xeqw / 10))
                            for i in results[-2:]
                        ]
                    ):  # Checking if last two runs gave a valid result
                        trials.append(xeqw * trials[-1] / results[-1])
                        logger.debug(
                            " Using linear approximation for new guess: {0:e}".format(
                                trials[-1]
                            )
                        )
                    else:
                        trials.append(Secant(trials, results, xeqw))
                        if trials[-1] < 0:
                            trials[-1] = xeqw * trials[-2] / results[-1]
                            logger.debug(
                                " Using linear approximation for new guess: {0:e}".format(
                                    trials[-1]
                                )
                            )
                        else:
                            logger.debug(
                                " Using secant method for new guess: {0:e}".format(
                                    trials[-1]
                                )
                            )
                if trials[-1] > 1.0:
                    logger.warning("Line Strength Insufficient")
                    finAbun.append(
                        {
                            "result": "error",
                            "message": "Line Strength Insufficient. Manual Examination suggested.",
                        }
                    )
                    break
            else:
                alleqw = CalcEqw(testLine)[1]
                wingpercent = ((alleqw - allzero) / results[-1] - 1) * 100
                finAbun.append(
                    {
                        "result": "success",
                        "relabun": trials[-1],
                        "logabun": math.log(trials[-1], 10) + getconf("LOGATREF"),
                        "wingpercent": wingpercent,
                    }
                )
                # finAbun.append('{0: >8.2e}  {1: >7.2f}   {2: >4.0f}%'.format(trials[-1], math.log(trials[-1],10) + getconf('LOGATREF'), wingpercent))
            logger.info(
                "Result: %s",
                (
                    finAbun[-1]["relabun"]
                    if finAbun[-1]["result"] == "success"
                    else finAbun[-1]["message"]
                ),
            )

        # Writing the output
        outputData = (
            {
                "model": args.model,
                "temperature": IS.TEMP,
                "logg": IS.LOGG,
                "version": version,
            },
            [],
        )
        if "unit55" in conf:
            outputData[0]["unit55"] = {
                param.upper(): getattr(IS, param.upper())
                for param in conf["unit55"]
                if hasattr(IS, param.upper())
            }
        for i, tl in enumerate(testLines):
            if type(tl) == str:
                outputData[1].append({"type": "comment", "value": tl.rstrip("\n")})
                continue
            outputData[1].append(
                {
                    "type": "line",
                    "target": tl[1],
                    "abundance": finAbun[i],
                    "lines": [
                        {
                            "wavelength": line.ALAM,
                            "ion": "{0: >2d}.{1:0>2d}".format(line.Z, line.Q),
                        }
                        for line in tl[0]
                    ],
                }
            )
        logger.debug("Writing Output")
        outputformatter[getconf("OUTFMT")](outputData, getconf("OUTFN"))
        logger.info("Total runs: %d", IS.runs)
except aeqwISError as e:
    if __debug__:
        raise


logger.info("Runtime: {0:.3f}".format(time() - startTime))
