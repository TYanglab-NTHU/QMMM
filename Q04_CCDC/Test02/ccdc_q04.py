# 13.08.2024 - Testing molecule from Ting Yee
# Read precisely one molecule and start to analyze it

# Load modules for Python coding
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import subprocess               # call the shell

# Load modules and prepare to read the CCS data base
from ccdc          import io                    # read crystal structures
from ccdc.search   import TextNumericSearch     # search by text or numbers
import ccdc.search

from ccdc.search import SubstructureSearch

from ccdc.molecule import Molecule, Atom, Bond  # build a molecule
csd_reader = io.EntryReader('CSD')

################################################################################
# Functions defined by me (might go into a module on their own)                #
################################################################################

# Fake printf and fprintf as used in C
def printf(format, *args):
  sys.stdout.write(format % args)
def fprintf(stream, format, *args):
  stream.write(format % args)

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

# load a xyz file
def load_xyz(name, ident):
  printf("Reading data from '%s'", name)
  # read all lines at once
  with open(name, 'r') as file:
    lines = file.readlines()
  NumAt=int(lines[0].strip())
  printf("  ...  %i atoms read\n", NumAt)
  # create an empty molecule
  mol = Molecule(identifier=ident)
  # add atoms line by line
  for line in lines[2:]:
    parts = line.split()
    at = Atom(parts[0], coordinates=(float(parts[1]), float(parts[2]), float(parts[3])))
    mol.add_atom(at)
  # try to bring some chemical meaning by brute force
  if bool(True):
    printf("  . Add missing H atoms\n")
    mol.add_hydrogens()
    printf("  . Identify rings\n")
    mol.assign_bond_types()
    mol.standardise_aromatic_bonds()
    mol.standardise_delocalised_bonds()
    # print(mol.bonds)
    
    print(mol.rings)
  return(mol)

################################################################################
# main code of the test                                                        #
################################################################################

# clean up old files
printf("Clean up old files\n")
if bool(False):
  printf("Clean up old files\n")
  printf("  Delete qm_at*xyz")
  for file in glob.glob("qm_at*xyz"):
    os.remove(file)
  printf("  Delete combi.xyz\n")
  os.remove("combi.xyz")
  
# read xyz file
if bool(False):
  mol = load_xyz("large.xyz", "test molecule from Ting Yee")
  write_xyz(mol, "back.xyz", "write back file")
else:
  if not os.path.isfile("./obt.mol2"):
    printf("Use OpenBabel to create obt.mol2\n")
    subprocess.run(["obabel", "./large.xyz", "-O", "./obt.mol2"])
  printf("read obt.mol2\n")
  mol_reader = io.MoleculeReader("./obt.mol2")
  mol = mol_reader[0]
  mol.assign_bond_types()
  mol.standardise_aromatic_bonds()
  mol.standardise_delocalised_bonds()

# A simple test for aromatic rings
if bool(False) and len(mol.rings)>0:
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

################################################################################
# identify functional groups                                                   #
################################################################################

# list of identfied functional groups to used later
func_groups=[]

printf("Search for functional groups\n")
# list of smart strings for the search of functional groups
# later rewrite this code to read the smarts form file
smarts_list=[]
smarts_list.append("C#N")         # nitrile / isonitile
# smarts_list.append("ccccccC#N")   # benzonitrile
# Cl atoms will not be added, because they link to the metal
smarts_list.append("cCl")         # Cl-C_ar (pos. test)
smarts_list.append("N#N")         # N2 gas  (neg. test)

# loop over list of smarts for functional groups
for smarts in smarts_list:
  # create a substructure form the smarts string
  substrcuct = ccdc.search.SMARTSSubstructure(smarts)
  # search for the substructure
  substructure_search = SubstructureSearch()
  substructure_search.add_substructure(substrcuct)
  matches = substructure_search.search(mol)
  # Check and report matches
  if matches:
    print(f'  Found {len(matches)} {smarts} group(s) in the molecule.')
    cnt=0
    for match in matches:
      printf("    %2i", cnt)
      fgroup=[]
      for at in match.match_atoms():
        str = "%s-%i" % (at.label, at.index)
        printf("  %-6s", str)
        fgroup.append(at)
      printf("\n")
      cnt+=1
      func_groups.append(fgroup)
  else:
    print(f'  No {smarts} group(s) found.')
print(f'  {len(smarts)} functional group(s) found.')
 
################################################################################
# the actual separation code                                                   #
################################################################################

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
  if len(ne.rings) != 0:
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

# the search for ligand rings ends here
# the search for ligand chains begins
# add functional groups

qm_cnt += 1
shell03=[]
for ne_at in shell01: # Loop over direct neighbors
  for fg in func_groups:
    if bool(False):
      printf("%s-%i", ne_at.label, ne_at.index)
      printf(" [")
      for z in fg:
        printf("%s-%i ", z.label, z.index)
      printf("]\n")
    if ne_at in fg:
      for at in fg:
        shell03.append(at)
shell03 = list(set(shell03))
# print(shell03)
fname = "qm_at%02i.xyz" % qm_cnt
list_xyz(shell03, fname, "Center, func. groups")

# add the neighborhood of the functional group atoms
# build a list of all atoms neighboring the functional groups

qm_cnt += 1
new_at = []
for fgat in shell03:
  for at in fgat.neighbours:
    if at not in shell02:
      new_at.append(at)
      if len(at.rings) != 0:
        if at.rings[0].is_aromatic:
          for z in at.rings[0].atoms:
            new_at.append(z)
if len(new_at) > 0:
  new_at=list(set(new_at))
  shell03 += new_at
  shell03=list(set(shell03))
  del z
  del new_at
fname = "qm_at%02i.xyz" % qm_cnt
list_xyz(shell03, fname, "Center, func. groups, neighbors")
# merge the results
qm += shell03
qm=list(set(qm))

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