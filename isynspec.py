# Isynspec.py
# Program to automate the calculation of abundance of an element given its
# equivalent width in a particular spectral line.
# K.Sriram
# Created: 20/04/2017

from subprocess import call

class InvalidInput(Exception):
    def __init__(self,line):
        self.msg = line
    def __str__(self):
        return self.msg

class INLIN:
    ALAM, Z, Q, GF, EXCL, QL, EXCU, QU, GS, GW, INEXT = 0., 1, 0, 0.,  0., 0., 0., 0.,  0., 0., 0.,
    def __init__(self, INLINstr):
        try:
            tokens = INLINstr.split()
            self.ALAM = float(tokens[0])
            ANUM = tokens[1].split('.')
            self.Z, self.Q = int(ANUM[0]), int(ANUM[1])
            self.GF = float(tokens[2])
            self.EXCL = float(tokens[3])
            self.QL = float(tokens[4])
            self.EXCU = float(tokens[5])
            self.QU = float(tokens[6])
            self.GS = float(tokens[7])
            self.GW = float(tokens[8])
            self.INEXT = float(tokens[9])
        except Exception as e:
            print e
            raise InvalidInput(INLINstr)
    def __str__(self):
        return '{0:.4f} {1:2d}.{2:0>2d}  {3:.3f}  {4:.3f} {5:.1f}  {6:.3f} {7:.1f}    {8:.2f}  {9:.2f}  {10:.2f} 0\n'.format(self.ALAM, self.Z, self.Q, self.GF, self.EXCL, self.QL, self.EXCU, self.QU, self.GS, self.GW, self.INEXT)
    def __cmp__(self,other):
        return self.ALAM - other.ALAM
        


class ISynspec:
    runs = 0
    # Here come some default values of all the parameters. See synspec guide to understand.
    modelFN = ''
    # fort.55
    IMODE, IDSTD, IPRIN = 1, 32, 0
    INMOD, INTRPL, ICHANG, ICHEMC = 1, 0, 0, 1
    IOPHLI = 0
    IFREQ, INLTE, ICONTL, INLIST, IFHE2 = 1, 1, 0, 0, 0
    IHYDPR, IHE1PR, IHE2PR = 0, 0, 0
    ALAM0, ALAM1, CUTOF0, CUTOFS, RELOP, SPACE = 2000.0, 5000.0, 15, 50, 0.000000000100, 0.01
    VTB = 13.0
    # fort.19
    LINELIST = []
    # fort.56
    ABUNDANCES = [] # To be stored as a list of tuples (int, float)
    # The template for fort.55
    temp55 = ( '{0:d} {1:d} {2:d}\n'
               '{3:d} {4:d} {5:d} {6:d}\n'
               '{7:d}\n'
               '{8:d} {9:d} {10:d} {11:d} {12:d}\n'
               '{13:d} {14:d} {15:d}\n'
               '{16:.1f} {17:.1f} {18:d} {19:d} {20:.12f} {21:f}\n'
               '{22:f}\n')
    def __init__(self):
        pass
    # Test if the Temperature and Gravity are same in both the model and concentration file
    def TestTG(self):
        with open('fort.8') as f:
            tokens = f.readline().split()
            self.TEMP = float(tokens[2])
            self.LOGG = float(tokens[3])
        with open(self.modelFN) as f:
            tokens = f.readline().split()
            TEMP = float(tokens[0])
            LOGG = float(tokens[1])
        return TEMP == self.TEMP and LOGG == self.LOGG

    # Methods to write to input files.
    def write55(self):
        with open('fort.55','w') as f:
            f.write(self.temp55.format(self.IMODE,self.IDSTD,self.IPRIN,self.INMOD,self.INTRPL,self.ICHANG,self.ICHEMC,self.IOPHLI,self.IFREQ,self.INLTE,self.ICONTL,self.INLIST,self.IFHE2,self.IHYDPR,self.IHE1PR,self.IHE2PR,self.ALAM0,self.ALAM1,self.CUTOF0,self.CUTOFS,self.RELOP,self.SPACE,self.VTB))
    def write19(self):
        with open('fort.19','w') as f:
            for line in self.LINELIST:
                f.write(str(line))
    def write56(self):
        with open('fort.56','w') as f:
            f.write('{0:d}\n'.format(len(self.ABUNDANCES)))
            for ABUN in self.ABUNDANCES:
                f.write('{0:d} {1:e}\n'.format(ABUN[0],ABUN[1]).replace('e',''))
    # Methods to read from output file
    def read16(self):
        self.EQW = []
        with open('fort.16') as f:
            for line in f:
                tokens = line.split()
                self.EQW.append( ( (float(tokens[0]),float(tokens[1])), float(tokens[2]) ) )
        return self.EQW
    # Running the SYNSPEC program. self.runs is a counter that keeps track of number of runs.
    def run(self):
        self.runs += 1
        with open(self.modelFN) as inp:
            out = open('/dev/null','w')
            p = call('./synbaba',stdin=inp,stdout=out)


