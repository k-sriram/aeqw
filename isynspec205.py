# Isynspec.py
# Program to automate the calculation of abundance of an element given its
# equivalent width in a particular spectral line.
# K.Sriram
# Created: 20/04/2017

from subprocess import call
import logging

logger = logging.getLogger('aeqw.iSynspec')

def fortfloat(x):
    if 'e' in x:
        return float(x)
    return float(x.replace('-','e-').replace('+','e+'))

class InvalidInput(Exception):
    def __init__(self,line):
        self.msg = line
    def __str__(self):
        return self.msg

class INLIN:
    ALAM, Z, Q, GF, EXCL, QL, EXCU, QU, GS, GW, INEXT = 0., 1, 0, 0.,  0., 0., 0., 0.,  0., 0., 0.,
    def __init__(self, INLINstr):
        try:
            line = ''.join([' '] * (10 - len(INLINstr.split()[0])) + [INLINstr.lstrip()])
            self.ALAM = float(line[0:10])
            self.Z, self.Q = int(line[11:13]), int(line[14:16])
            self.GF = float(line[16:23])
            self.EXCL = float(line[23:35])
            self.QL = float(line[35:39])
            self.EXCU = float(line[39:51])
            self.QU = float(line[51:55])
            self.AGAM = float(line[55:63])
            self.GS = float(line[63:70])
            self.GW = float(line[70:77])
            self.remainder = line[79:]
        except Exception as e:
            logger.error("Error parsing spectral line: ")
            logger.error(INLINstr)
            logger.error(e)
            raise InvalidInput(INLINstr)
    def __str__(self):
        return '{0:>10.4f}{1:>3d}.{2:0>2d}{3:>7.3f}{4:>12.3f}{5:>4.1f}{6:>12.3f}{7:>4.1f} {8:>7.2f}{9:>7.2f}{10:>7.2f} 0'.format(self.ALAM, self.Z, self.Q, self.GF, self.EXCL, self.QL, self.EXCU, self.QU, self.AGAM, self.GS, self.GW)
        


class ISynspec:
    runs = 0
    # Here come some default values of all the parameters. See synspec guide to understand.
    modelFN = ''
    # fort.55
    IMODE, IDSTD, IPRIN = 1, 32, 0
    INMOD, INTRPL, ICHANG, ICHEMC = 1, 0, 0, 1
    IOPHLI, NUNALP, NUNBET, NUNGAM, NUNBAL = 0, 0, 0, 0, 0
    IFREQ, INLTE, ICONTL, INLIST, IFHE2 = 1, 1, 0, 0, 0
    IHYDPR, IHE1PR, IHE2PR = 0, 0, 0
    CUTOF0, RELOP, SPACE = 15, 0.000000000100, 0.01
    ALAM0, ALAM1, CUTOFS = 2000.0, 5000.0, 50
    NMLIST = '0 0'
    VTB = 13.0
    # fort.19
    LINELIST = []
    # fort.56
    ABUNDANCES = [] # To be stored as a list of tuples (int, float)
    # The template for fort.55
    temp55 = ( '{0:d} {1:d} {2:d}\n'
               '{3:d} {4:d} {5:d} {6:d}\n'
               '{7:d} {8:d} {9:d} {10:d} {11:d}\n'
               '{12:d} {13:d} {14:d} {15:d} {16:d}\n'
               '{17:d} {18:d} {19:d}\n'
               '{20:.1f} {21:.1f} {22:d} {23:d} {24:.1e} {25:f}\n'
               '{26}\n'
               '{27:f}\n')
    def __init__(self):
        self.read55()
    # Test if the Temperature and Gravity are same in both the model and concentration file
    # An error was found in this method of reading fort.8. It seems not all fort.8 files have this information.
    def TestTG(self):
        logger.debug(" Checking to see if T and log g are consistent between files.")
        with open(self.modelFN+'.5') as f:
            tokens = f.readline().split()
            self.TEMP = float(tokens[0])
            self.LOGG = float(tokens[1])
        try:
            with open(self.modelFN+'.7') as f:
                tokens = f.readline().split()
                TEMP = float(tokens[2])
                LOGG = float(tokens[3])
        except IndexError:
            logger.warning("Unable to read temperature and gravity information in fort.8")
            return True # Reurning True although there is an error since there is no mismatch. Only absence of data.
        else:
            return TEMP == self.TEMP and LOGG == self.LOGG

    # Methods to write to input files.
    def write55(self):
        logger.debug("   Writing to fort.55")
        with open('fort.55','w') as f:
            f.write(self.temp55.format(self.IMODE,self.IDSTD,self.IPRIN,self.INMOD,self.INTRPL,self.ICHANG,self.ICHEMC,self.IOPHLI,self.NUNALP,self.NUNBET,self.NUNGAM,self.NUNBAL,self.IFREQ,self.INLTE,self.ICONTL,self.INLIST,self.IFHE2,self.IHYDPR,self.IHE1PR,self.IHE2PR,self.ALAM0,self.ALAM1,self.CUTOF0,self.CUTOFS,self.RELOP,self.SPACE,self.NMLIST,self.VTB).replace('e',''))
    def write19(self):
        logger.debug("   Writing to fort.19")
        with open('fort.19','w') as f:
            for line in self.LINELIST:
                f.write(str(line))
                f.write('\n')
    def write56(self):
        logger.debug("   Writing to fort.56")
        with open('fort.56','w') as f:
            f.write('{0:d}\n'.format(len(self.ABUNDANCES)))
            for ABUN in self.ABUNDANCES:
                # Checking for unusual abundances
                if ABUN[1] > 1.0 or ABUN[1] < 0.0:
                    logger.warning("Unusual value of abundance is being input: {0:d} {1:e}".format(ABUN[0],ABUN[1]))
                f.write('{0:d} {1:e}\n'.format(ABUN[0],ABUN[1]).replace('e',''))
    # Methods to read from output file.
    def read16(self):
        logger.debug("   Reading from fort.16")
        self.EQW = []
        with open('fort.16') as f:
            for line in f:
                tokens = line.split()
                self.EQW.append( ( (float(tokens[0]),float(tokens[1])), float(tokens[2]) ) )
        return self.EQW
    # Reading fort.55 for values of all parameter.
    def read55(self):
        logger.debug("    Reading from fort.55")
        with open('fort.55') as f:
            # Line 1
            self.IMODE, self.IDSTD, self.IPRIN = [int(i) for i in f.readline().split()]
            # Line 2
            self.INMOD, self.INTRPL, self.ICHANG, self.ICHEMC = [int(i) for i in f.readline().split()]
            # Line 3
            self.IOPHLI, self.NUNALP, self.NUNBET, self.NUNGAM, self.NUNBAL = [int(i) for i in f.readline().split()]
            # Line 4
            self.IFREQ, self.INLTE, self.ICONTL, self.INLIST, self.IFHE2 = [int(i) for i in f.readline().split()]
            # Line 5
            self.IHYDPR, self.IHE1PR, self.IHE2PR = [int(i) for i in f.readline().split()]
            # Line 6
            tokens = f.readline().split()
            self.CUTOF0, self.RELOP, self.SPACE = int(tokens[2]), fortfloat(tokens[4]), float(tokens[5])
            # Line 7
            self.NMLIST = f.readline().strip()
            # Line 8
            self.VTB = float(f.readline().strip())
    # Reading fort.56, this can be used to set the initial abundances or checking.
    def read56(self):
        logger.debug("    Reading from fort.56")
        with open('fort.56') as f:
            self. ABUNDANCES = []
            lines = int(f.readline().strip())
            for i in range(lines):
                line = f.readline()
                tokens = line.split()
                self.ABUNDANCES.append((int(tokens[0]),fortfloat(tokens[1])))


    # Running the SYNSPEC program. self.runs is a counter that keeps track of number of runs.
    def run(self):
        logger.debug("   Running SYNSPEC: RSynspec {0}".format(self.modelFN))
        call(['rm','-f','fort.16']) # To avoid reading previous data in case of SYNSPEC not running.
        self.runs += 1
        with open('/dev/null','r+') as nullf:
            call(['RSynspec', self.modelFN],stdout=nullf)


