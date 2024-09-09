# 08.08.2024 - Getting started
# Read precisely one molecule and start to analyze it
# The identifier is CIMHIT
# https://www.ccdc.cam.ac.uk/structures/Search?Doi=10.1038%2Fs41586-023-05821-2&DatabaseToSearch=Published
# I postpone the problem of disordered crystal structures and continue with my ad hoc solution

# Load modules for Python coding
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system

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
  printf("Unable to retrieve CCDB entry for <%s>\n", ident)
  exit()
else:
  printf("Retrieving CCDB entry for <%s>\n", ident)
try:
  mol = csd_reader.molecule(ident)
except:
  printf("Unable to retrieve molecular geometry for <%s>\n", ident)
  exit()
else:
  printf("Retrieving molecular geometry for <%s>\n", ident)
  
# write a mol2 file to check the integrety of the geo data
cnt = 1
fname = "entry%02i.mol2" % cnt
io.MoleculeWriter(fname).write(mol)
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

# pick the metal center
cent = []
qm = []
mm =[]
for at in mol.atoms:
  printf("%-4s  ", at.label)
  if at.is_metal:
    printf("M  ")
    if at.atomic_number >= 19:
      printf("H  ")
      cent.append(at)
      qm.append(at)
    else:
      printf("-  ")
  else:
    printf("-  -  ")
  printf("\n")
if len(cent) < 1:
  printf("No metal center found\n")
  exit()
else:
  printf("%i metal center found\n", len(cent))
printf("\n")

# the center and its neighbours
printf("The metal center and its neighbours\n")
neighbours = []
for at in cent:
  printf("%-4s  ", at.label)
  for na in at.neighbours:
    printf("%-4s  ", na.label)
    neighbours.append(na)
    qm.append(na)
  printf("\n")
qm = list(set(qm))
printf("\n")
 
# check if the neighbours are part of ring
printf("Are the neighbours part of a ring?\n")
for na in neighbours:
  printf("%-4s  ", na.label)
  if na.is_cyclic:
    printf("Y  ")
    for ra in na.rings[0].atoms:
      printf("%-4s  ", ra.label)
      qm.append(ra)
  else:
    printf("n  ")
  printf("\n")
qm = list(set(qm))

# write qm atoms to a file
out = open("qm.xyz", "w")
fprintf(out, "%s\n", len(qm))
fprintf(out, "%s\n", "Quantum atoms")
for at in qm:
  if hasattr(at.coordinates, "x"):
    fprintf(out, "%-2s  ",   at.atomic_symbol)
    fprintf(out, "%10.6f  ", at.coordinates.x)
    fprintf(out, "%10.6f  ", at.coordinates.y)
    fprintf(out, "%10.6f\n", at.coordinates.z)
out.close()




exit()