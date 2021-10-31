#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# autoeqw.py
# Program to automate the calculation of abundance of an element given its
# equivalent width in a particular spectral line.
# K.Sriram
# Created: 20/04/2017


import sys
from time import time
from datetime import datetime
import logging
import logging.handlers
import math
from configparser import ConfigParser, ExtendedInterpolation
from argparse import ArgumentParser
import json
from aeqw.isynspec import ISynspec, INLIN
from aeqw import __version__

CONFFN = "aeqw.conf"


def parse_cmd(argv=None):
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
        "-c",
        action="append",
        help="Extra configuration options files. Can be repeated.",
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
    argparser.add_argument(
        "--logatref", type=float, help="Abundance of reference element."
    )
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
    argparser.add_argument(
        "-V", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    args = argparser.parse_args(argv)
    return args


argconf = (
    "infn",
    "outfn",
    "extralogfn",
    "initabun",
    "nullabun",
    "logatref",
    "broad",
    "range",
    "epsilon",
    "outfmt",
)
argconfbool = ("sep19",)


class Config(ConfigParser):
    confgetter = {
        "str": "get",
        "float": "getfloat",
        "int": "getint",
        "bool": "getboolean",
    }

    def __init__(self, defaultconf):
        super().__init__(interpolation=ExtendedInterpolation())
        self["DEFAULT"] = {
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
        self["TYPES"] = {
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
        self["aeqw"] = {}
        self.sec = "aeqw"
        self.readconf(defaultconf)

    def readconf(self, conffn, warnifnotfound=False):
        if self.read(conffn) == []:
            if warnifnotfound:
                logger.warn(f" Configuration file {conffn} not found.")
            else:
                logger.info(f" Configuration information can be added to file {conffn}")
            return
        logger.info(f" Reading Configuration file: {conffn}")

    def addconfs(self, confs):
        if confs is not None:
            for c in confs:
                self.readconf(c, True)

    def getconf(self, param, sec=None):
        if sec is None:
            sec = self.sec
        if param in self[sec]:
            try:
                return getattr(self[sec], self.confgetter[self["TYPES"][param]])(param)
            except ValueError:
                logger.error(f"Invalid type for parameter {param} in configuration")
                raise
        else:
            logger.error(f"Parameter {param} not found")
            raise ValueError

    def add_args(self, args, argconf, argconfbool):
        for c in argconf:
            if getattr(args, c) is not None:
                self["aeqw"][c.upper()] = getattr(args, c)
        for c in argconfbool:
            if getattr(args, c) == True:
                self["aeqw"][c.upper()] = True


def Overlap(bin, box):
    l = max(bin[0], box[0])
    r = min(bin[1], box[1])
    if l > r:
        return 0.0
    return math.sqrt((r - l) / (bin[1] - bin[0]))


def Secant(x, f, y, epsilon):
    if abs(f[-1] - f[-2]) < epsilon:
        return -1
    return (x[-2] * (f[-1] - y) - x[-1] * (f[-2] - y)) / (f[-1] - f[-2])


# Output formatters


def outputtxt(outputData, outfn):
    with open(outfn, "w") as f:
        f.write("{model:s} {temperature:.2f} {logg:.2f}\n".format_map(outputData[0]))
        if "unit55" in outputData[0]:
            f.write(
                "".join(
                    [f"{key} = {value}\n" for key, value in outputData[0]["unit55"]]
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
                f.write(f" {row['target']:8.2f}  {abuntxt}")
        f.write("\n")


def outputjson(outputData, outfn):
    txt = json.dumps(outputData)
    with open(outfn, "w") as f:
        f.write(txt)


outputformatter = {"txt": outputtxt, "json": outputjson}


def aeqw(conf, model, outputformatter):
    with ISynspec(model) as synspec_interface:

        if "unit55" in conf:
            for param in conf["unit55"]:
                if hasattr(synspec_interface, param.upper()):
                    setattr(
                        synspec_interface,
                        param.upper(),
                        type(getattr(synspec_interface, param.upper()))(
                            conf["unit55"][param]
                        ),
                    )
                    logger.info(
                        f"Setting unit 55 parameter {param.upper()} to {getattr(synspec_interface, param.upper())}"
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
        with open(conf.getconf("INFN")) as f:
            for line in f:
                logger.debug(f" Processing line: {line.strip()}")
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
                    logger.error(f"Error while processing line {lineNo:d}\n{line}\n")
                    raise

        allLines.sort(key=lambda x: x.ALAM)
        synspec_interface.LINELIST = allLines
        synspec_interface.write19()
        # Determining the bounds of the synthetic spectrum ALAM0 and ALAM1. The multiplication by ten is for conversion from nm to A. Also writing to 19 and 55

        def InitParam(testLine, allLines):
            synspec_interface.ALAM0 = min(
                [line.ALAM for line in testLine]
            ) * 10 - conf.getconf("RANGE")
            synspec_interface.ALAM1 = max(
                [line.ALAM for line in testLine]
            ) * 10 + conf.getconf("RANGE")
            logger.debug(
                f" InitParam: Setting range of synthetic spectrum: ({synspec_interface.ALAM0:.1f}, {synspec_interface.ALAM1:.1f})"
            )

            synspec_interface.write55()

            if conf.getconf("SEP19") == True:
                synspec_interface.LINELIST = testLine
                synspec_interface.write19()

        # Calculate the Equivalent width of a particular line
        def CalcEqw(testLine):
            if len(synspec_interface.EQW) < 2:
                logger.warning("  CalcEqw: SYNSPEC did not generate output in fort.16")
                return None, 0
            box = (
                min([line.ALAM for line in testLine]) * 10 - conf.getconf("BROAD"),
                max([line.ALAM for line in testLine]) * 10 + conf.getconf("BROAD"),
            )
            logger.debug(
                f"  CalcEqw: Calculating Equivalent width; including bins in {box}."
            )
            total = 0
            alltotal = 0
            for bin in synspec_interface.EQW:
                total += bin[1] * Overlap(bin[0], box)
                alltotal += bin[1]
            logger.debug(f"  CalcEqw: eqw = {total:f}, alleqw = {alltotal:f}")
            return total, alltotal

        # Set the abundance and run SYNSPEC and read the output
        def Run(abundances):
            logger.debug(f"  Setting abundance: {abundances}")
            synspec_interface.ABUNDANCES = abundances
            synspec_interface.write56()
            synspec_interface.run()
            synspec_interface.read16()

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
            abun = conf.getconf("NULLABUN")
            Run([(Z, abun)])
            zero, allzero = CalcEqw(testLine)
            while zero is None:  # If the program didn't compute the bins
                if synspec_interface.RELOP > 1e-12:
                    synspec_interface.RELOP /= 10
                    logger.debug(
                        f" > Setting RELOP parameter to {synspec_interface.RELOP:.1e}"
                    )
                    synspec_interface.write55()
                Run([(Z, abun)])
                zero, allzero = CalcEqw(testLine)

            logger.debug(f" > Zero = {zero:f}, allZero = {allzero:f}")

            # Finding the abundance that gives reasonable eqw
            trials = [synspec_interface.INITABUNZWISE.get(Z, conf.getconf("INITABUN"))]
            results = []
            while not results or abs(results[-1] - xeqw) > conf.getconf("EPSILON"):
                logger.debug(
                    f" Running for abundance: {trials[-1]:e}, target width: {xeqw:f}"
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
                        f"  Guess = {trials[-1]:e}, Result = {results[-1]:f}, Target = {xeqw:f}, Diff = {xeqw - results[-1]:f}, Epsilon = {conf.getconf('EPSILON'):f}"
                    )
                    if len(results) < 2 or all(
                        [
                            not (i is None or (i >= 0 and i < xeqw / 10))
                            for i in results[-2:]
                        ]
                    ):  # Checking if last two runs gave a valid result
                        trials.append(xeqw * trials[-1] / results[-1])
                        logger.debug(
                            f" Using linear approximation for new guess: {trials[-1]:e}"
                        )
                    else:
                        trials.append(
                            Secant(trials, results, xeqw, conf.getconf("EPSILON"))
                        )
                        if trials[-1] < 0:
                            trials[-1] = xeqw * trials[-2] / results[-1]
                            logger.debug(
                                f" Using linear approximation for new guess: {trials[-1]:e}"
                            )
                        else:
                            logger.debug(
                                f" Using secant method for new guess: {trials[-1]:e}"
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
                        "logabun": math.log(trials[-1], 10) + conf.getconf("LOGATREF"),
                        "wingpercent": wingpercent,
                    }
                )
                # finAbun.append(f"{trials[-1]: >8.2e}  {math.log(trials[-1],10) + conf.getconf('LOGATREF'): >7.2f}   {wingpercent: >4.0f}%")
            logger.info(
                f"Result: {finAbun[-1]['relabun'] if finAbun[-1]['result'] == 'success' else finAbun[-1]['message']}"
            )

        # Writing the output
        outputData = (
            {
                "model": model,
                "temperature": synspec_interface.TEMP,
                "logg": synspec_interface.LOGG,
                "version": __version__,
            },
            [],
        )
        if "unit55" in conf:
            outputData[0]["unit55"] = {
                param.upper(): getattr(synspec_interface, param.upper())
                for param in conf["unit55"]
                if hasattr(synspec_interface, param.upper())
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
                        {"wavelength": line.ALAM, "ion": f"{line.Z: >2d}.{line.Q:0>2d}"}
                        for line in tl[0]
                    ],
                }
            )
        logger.debug("Writing Output")
        outputformatter(outputData, conf.getconf("OUTFN"))
        logger.info("Total runs: %d", synspec_interface.runs)


def init_logger():
    global logger
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


def add_extralog(args, conf):
    extralog = (
        args.extralogfn
        if args.extralogfn is not None
        else (conf.getconf("EXTRALOGFN") if conf.getconf("EXTRALOGFN") != "" else None)
    )

    if extralog is not None:
        xtrafilelog = logging.FileHandler(extralog)
        xtrafilelog.setLevel(logging.DEBUG)
        xtrafilelog.setFormatter(
            logging.Formatter("%(asctime)s - %(name)-15s %(levelname)-8s: %(message)s")
        )
        logger.addHandler(xtrafilelog)


def main(argv=None):
    startTime = time()

    args = parse_cmd(argv)
    model = args.model

    init_logger()
    logger.info(
        f"Running program: Automatic Equation width solver Version: {__version__}"
    )
    logger.info(str(datetime.now()))

    conf = Config(CONFFN)

    add_extralog(args, conf)

    conf.add_args(args, argconf, argconfbool)

    logger.info(f"Model: {model}")
    logger.debug("Initialization")

    logger.debug(" Parameters:")
    for param in conf["aeqw"].keys():
        logger.debug(f"  {param} : {conf.getconf(param)}")

    aeqw(conf, model, outputformatter[conf.getconf("OUTFMT")])

    logger.info(f"Runtime: {time() - startTime:.3f}")


if __name__ == "__main__":
    main()
