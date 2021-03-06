'''
Container for chemical data
Created on Jan 30, 2011

@author Jiahao Chen
'''

## A periodic table with element names
Elements = [x.upper() for x in \
           ['H' , 'He',
            'Li', 'Be', 'B' , 'C' , 'N', 'O' , 'F' , 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P', 'S' , 'Cl', 'Ar',
            'K' , 'Ca',
                  'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                       'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
            'Rb', 'Sr',
                  'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                        'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe',
            'Cs', 'Ba',
                  'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
                  'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
                  'Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
                        'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
            'Fr', 'Ra',
                  'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu', 'Am',
                  'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
                  'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn'
                       ]
           ]

#The working assumption is that everything in the code is handled in atomic
#units until it's time to print something output or interface
#The convention is to mimic in the code how it would be written,
#where possible, e.g.
#1.2 * eV, 0.13 * Hartree, etc.

Hartree = 1.0
eV = 0.03674932534 #converts eV to Hartree
Angstrom = 1.889725989
Debye = 1/2.5417462
HartreeToNm = 1239.8*eV

## If run as main, prints out all the elements
if __name__ == '__main__':
    for n, e in enumerate(Elements):
        print n + 1, e
