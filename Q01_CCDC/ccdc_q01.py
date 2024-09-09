# 08.08.2024 - Getting started
# Read precisely one molecule and start to analyze it
# The identifier is CIMHIT
# https://www.ccdc.cam.ac.uk/structures/Search?Doi=10.1038%2Fs41586-023-05821-2&DatabaseToSearch=Published

# Load modules for Python coding
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system

# writing the xyz file from molecule
def write_xyz(mol, name, remark):
  out = open(name, "w")
  fprintf(out, "%s\n", len(mol.atoms))
  fprintf(out, "%s\n", remark)
  for n in range(len(mol.atoms)):
    fprintf(out, "%-2s  ", mol.atoms[n].atomic_symbol)
    fprintf(out, "%10.6f  ", mol.atoms[n].coordinates.x)
    fprintf(out, "%10.6f  ", mol.atoms[n].coordinates.y)
    fprintf(out, "%10.6f\n", mol.atoms[n].coordinates.z)
  out.close()
  return()

# writing the xyz file from a list of atoms
def list_xyz(alist, name, remark):
  out = open(name, "w")
  fprintf(out, "%s\n", len(alist))
  fprintf(out, "%s\n", remark)
  for n in range(len(alist)):
    fprintf(out, "%-2s  ", alist[n].atomic_symbol)
    fprintf(out, "%10.6f  ", alist[n].coordinates.x)
    fprintf(out, "%10.6f  ", alist[n].coordinates.y)
    fprintf(out, "%10.6f\n", alist[n].coordinates.z)
  out.close()
  return()

# clean up old files
for file in glob.glob("*mol2"):
  os.remove(file)
for file in glob.glob("*xyz"):
  os.remove(file)

# Load modules and prepare to read the CCS data base
from ccdc          import io                    # read crystal structures
from ccdc.search   import TextNumericSearch     # search by text or numbers
from ccdc.molecule import Molecule, Atom, Bond  # Build a molecule
csd_reader = io.EntryReader('CSD')

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
  
# write a mol2 file to check the integrety of the geo data
cnt = 1
fname = "start%02i.mol2" % cnt
io.MoleculeWriter(fname).write(mol)

# A simple test for aromatic rings
if bool(False):
  printf("xx  ")
  printf("conj   ")
  printf("arom   ")
  printf("ring atoms\n")
  rz=1
  for ri in mol.rings:
    printf("%2i  ", rz)
    printf("%-5s  ", ri.is_fully_conjugated)
    printf("%-5s", ri.is_aromatic)
    if ri.is_aromatic:
      printf("  %s\n", ri)
    else:
      printf("\n")
    rz += 1
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

# update structural information
if bool(True):
  printf("  . Assign unknown bonds\n")
  mol.assign_bond_types(which='unknown')
  printf("  . Standarize aromatic bonds\n")
  mol.standardise_aromatic_bonds()
  printf("  . Standarize delocalized bonds\n")
  mol.standardise_delocalised_bonds()
  printf("  . Add missing H atoms\n")
  mol.add_hydrogens()
    

# write the molecule to disc
cnt = 2 # increas counter to for the next step
fname = "start%02i.mol2" % cnt
io.MoleculeWriter(fname).write(mol)
fname = "start%02i.xyz" % cnt
write_xyz(mol, fname, "Cleaned up geometry")

# initialize lists
qm = [] # list with QM atoms
mm = [] # list with MM atoms
cent = [] # metal center
shell01 = [] # direct neighbours
shell02 = [] # rings to the direct neighbours

# find metal centers
qm_cnt = 1
for at in mol.atoms:
  if at.is_metal and at.atomic_number >= 19:
    cent.append(at)
    qm.append(at)
cent=list(set(cent))
qm=list(set(qm))
fname = "qm_at%02i.xyz" % qm_cnt
list_xyz(qm, fname, "Metal center")

# find the direct neighbours to the metal center
qm_cnt += 1
for mc in cent:
  for at in mc.neighbours:
    shell01.append(at)
    qm.append(at)
shell01=list(set(shell01))
qm=list(set(qm))
fname = "qm_at%02i.xyz" % qm_cnt
list_xyz(qm, fname, "Center, direct neighbours")

# add the atoms of rings containing direct neighbours
qm_cnt += 1
shell02 = shell01.copy()
for ne in shell01:
  for at in ne.rings[0].atoms: # focus on the smallest ring
    shell02.append(at)
    qm.append(at)
shell02=list(set(shell02))
qm=list(set(qm))
fname = "qm_at%02i.xyz" % qm_cnt
list_xyz(qm, fname, "Center, direct neighbour rings")

# are atoms in shell02 part of an aromatic system which needs to be added to shell02
qm_cnt += 1
new_at = []
for sa in shell02:
  for ri in sa.rings:
    if ri.is_aromatic:
      for at in ri.atoms:
        new_at.append(at)
shell02 += new_at
shell02=list(set(shell02))
qm += new_at
qm=list(set(qm))
fname = "qm_at%02i.xyz" % qm_cnt
list_xyz(qm, fname, "Center, aromatic rings")

# add H atoms connected to QM atoms
qm_cnt += 1
new_at = []
for at in qm:
  for na in at.neighbours:
    if na.atomic_number == 1:
      new_at.append(na)
new_at=list(set(new_at))
qm += new_at
qm=list(set(qm))
fname = "qm_at%02i.xyz" % qm_cnt
list_xyz(qm, fname, "neighbouring H atoms added")

# any atom not QM is MM
for at in mol.atoms:
  if at in qm:
    pass
  else:
    mm.append(at)
mm=list(set(mm))
fname = "mm_at00.xyz"
list_xyz(mm, fname, "MM atoms")

# combine both sets for testing
combi=[]
printf("Combined output:\n")
printf("  QM atoms: %i\n", len(qm))
printf("  MM atoms: %i\n", len(mm))
combi = qm+mm
fname = "combi.xyz"
comment = "QM: %i atoms, MM: %i atoms" % (len(qm), len(mm))
list_xyz(combi, fname, comment)

exit()