# 08.08.2024 - Getting started
# Read include funcionality to handle disorder
# The identifier is CIMHIT
# https://www.ccdc.cam.ac.uk/structures/Search?Doi=10.1038%2Fs41586-023-05821-2&DatabaseToSearch=Published

# Load modules for Python coding
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import numpy as np

# Load modules and prepare to read the CCS data base
from ccdc import io                       # read crystal structures
csd_reader = io.EntryReader('CSD')
from ccdc.search import TextNumericSearch # search by text or numbers

# Fake printf and fprintf as used in C
def printf(format, *args):
  sys.stdout.write(format % args)
def fprintf(stream, format, *args):
  stream.write(format % args)

# clean up old files
for file in glob.glob("entry*mol2"):
  os.remove(file)
for file in glob.glob("entry*xyz"):
  os.remove(file)
  
# Try to retrieve the structural infiotmation directly by the identifier
ident="CIMHIT"

# build a structure search by identifier
printf("Trying to read <%s>\n", ident)
try:
  ent = csd_reader.entry(ident)
except:
  printf("  Unable to retrieve CCDB entry for <%s>\n", ident)
  exit()
else:
  printf("  Retrieving CCDB entry for <%s>\n", ident)
try:
  mol = csd_reader.molecule(ident)
except:
  printf("  Unable to retrieve molecular geometry for <%s>\n", ident)
  exit()
else:
  printf("  Retrieving molecular geometry for <%s>\n", ident)
try:
  cryst = csd_reader.crystal(ident)
except:
  printf("  Unable to retrieve crystal structure for <%s>\n", ident)
  exit()
else:
  printf("  Retrieving crystal structure for <%s>\n", ident)

  
# write a mol2 file to check the integrety of the geo data
cnt = 1
fname = "entry%02i.mol2" % cnt
io.MoleculeWriter(fname).write(mol)

# the nitrile shows disorder, which needs to be fixed
printf("Checking crystal structure for disorder\n")
if cryst.has_disorder:
  printf("  Disorder detected -> Trying to clean up\n")
  molA = mol.copy()
  # removing all atoms with a label ending with a capital 'A' didn't help
  # I remove all A atoms closer than 0.7 A to an N atom
  # Differentiating by the label didn't help
  molN = mol.copy()
  for aA in molN.atoms:
    for aB in molN.atoms:
      if aA != aB:
        dist= np.linalg.norm(np.array(aA.coordinates) - np.array(aB.coordinates))
        if dist < 0.75:
          printf("%-4s  %-4s  %7.4f\n", aA.label, aB.label, dist)
          molN.remove_atom(aB)
  fname = "entry%02iN.mol2" % cnt
  io.MoleculeWriter(fname).write(molN)
else:
  printf("  No disorder\n")

exit()

# the geometry of the nitrile is definetly brocken, too many atoms

# in this specian case, I trry to fix the problem manually
dlabels  = [
  'C52', 'C53', 'C54', 'C55', 'C56', 'C57', 'H32', 'H33', 'H34', 'H35', 'H36', 'C51', 'N5',   # 2nd solvent molecule
  'N4', 'C44',                                                                                # CN mush
  'C45', 'C45A', 'C46', 'C46A', 'C47', 'C47A', 'C48','C48A', 'C49', 'C49A', 'C50', 'C50A',    # C6H5 mush
  'H27', 'H27A', 'H28', 'H28A', 'H29', 'H29A', 'H30', 'H30A', 'H31', 'H31A',                  # C6H5 mush
  'P1', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6',                                                   # PF6 anion
  ]
# delete the strange atoms
for at in mol.atoms:
  if at.label in dlabels:
    mol.remove_atom(at)
# write the molecule to disc
cnt = 2 # increas counter to for the next step
fname = "entry%02i.mol2" % cnt
io.MoleculeWriter(fname).write(mol)

# writing the final xyz file
fname = "entry%02i.xyz" % cnt
out = open(fname, "w")
fprintf(out, "%s\n", len(mol.atoms))
fprintf(out, "%s\n", ent.chemical_name)
for n in range(len(mol.atoms)):
  fprintf(out, "%-2s  ", mol.atoms[n].atomic_symbol)
  fprintf(out, "%10.6f  ", mol.atoms[n].coordinates.x)
  fprintf(out, "%10.6f  ", mol.atoms[n].coordinates.y)
  fprintf(out, "%10.6f\n", mol.atoms[n].coordinates.z)
out.close()

exit()