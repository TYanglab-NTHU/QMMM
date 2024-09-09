# 17.08.2024
#   - Testing bimetallic molecule from Ting Yee
#   - Read precisely one molecule and start to analyze it
# 18.08.2034
#   - Add a serach for briges between first layer rings
#   - Start reorganizing the code to preserve intermediate finds
# 19.08.2024
#   - Add aliphatic bridges between atoms in thelayer
# 20.08.2024
#   - Look for two atom links between rings
#   - Organize the shells in one list
#     * shell00 vontains the metal center for compability reasons

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
if bool(True):
  printf("Clean up old files\n")
  printf("  Delete qm_at*xyz\n")
  for file in glob.glob("qm_at*xyz"):
    os.remove(file)
  printf("  Delete mm_at01.xyz\n")
  if os.path.exists("./mm_at00.xyz"): os.remove("./mm_at00.xyz")
  printf("  Delete combi.xyz\n")
  if os.path.exists("./combi.xyz"): os.remove("./combi.xyz")
else:
  printf("Keep old files\n")

# read xyz file
if bool(False):
  mol = load_xyz("bimetal.xyz", "test Co/Sn molecule from Ting Yee")
  write_xyz(mol, "back.xyz", "write back file")
else:
  if not os.path.isfile("./obt.mol2"):
    printf("Use OpenBabel to create obt.mol2\n")
    subprocess.run(["obabel", "./bimetal.xyz", "-O", "./obt.mol2"])
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
# function InMolBonds to test if two atoms are joined by registered bond       #
# input  MyMol     molecule to be serach                                       #
#        At1       index of the first atom                                     #
#        At2       index of the second atom                                    #
# output -1        the bond is not MyMol.bonds                                 #
#        integer   index number in the bond in the list                        #
################################################################################

def InMolBonds(MyMol, At1, At2):
  SearchResult = -1  # default value
  for mb in MyMol.bonds:
    if bool(False):
      printf("  %3i", MyMol.bonds.index(mb))
      MyLabel = "%s-%i" % (mb.atoms[0].atomic_symbol, mb.atoms[0].index)
      printf("  %6s", MyLabel) 
      MyLabel = "%s-%i" % (mb.atoms[1].atomic_symbol, mb.atoms[1].index)
      printf("  %6s\n", MyLabel)
    if mb.atoms[0].index in [At1, At2] and mb.atoms[1].index in [At1, At2]:
      SearchResult = MyMol.bonds.index(mb)
  return(SearchResult)

################################################################################
# function FindSmarts to use SMARTS to parse the molecule                      #
# input  MySmarts  string to hold the SMARTS to serach for                     #
#        MyMol     molecule to searched                                        #
# output MyAtoms   a list of lists holding the identfied atoms                 #
################################################################################

def FindSmarts(MySmarts, MyMol):
  # inilize variables
  MyAtoms = []
  
  # create a substructure form the smarts string
  substrcuct = ccdc.search.SMARTSSubstructure(smarts)
  # search for the substructure
  substructure_search = SubstructureSearch()
  substructure_search.add_substructure(substrcuct)
  # serach for the substructure
  Matches = substructure_search.search(mol)
  
  # evalualte possible matches and build atom list
  if Matches:
    for MyMatch in Matches:
      AtomMatchList = []
      for at in MyMatch.match_atoms():
        AtomMatchList.append(at)
      MyAtoms.append(AtomMatchList)
      
  # return list with the results
  return(MyAtoms)

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
smarts_list.append("[#6R]A[#6R]")         # CH3 group (neg. test)

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
print(f'  {len(func_groups)} functional group(s) found.')
 
################################################################################
# the actual separation code                                                   #
################################################################################

# functions to help sorting lists
# sort a list of atoms by the index of the atom
def my_sort_by_index(ml):
  def my_key(ma):
    return(ma.index)
  ml.sort(key=my_key)
  return(ml)
  
# clean up my list and sort it by atom index
def my_clean_list(ml):
  ml=list(set(ml))     # remove dublicates
  ml=my_sort_by_index(ml) # sort the atom list by atom index
  return(ml)
  
# print a atom list 
def my_atom_list(ml, indent=""):
  for at in ml:
    printf("%s%3i  ",  indent, at.index)
    printf("%2s  ",    at.atomic_symbol)
    printf("%10.6f  ", at.coordinates.x)
    printf("%10.6f  ", at.coordinates.y)
    printf("%10.6f\n", at.coordinates.z)

    
# lists with intermediate results
# the Z indicates lists with intermediate results
Zneighbors  = [] # direct neighbours to the metal centers
Zcloserings = [] # rings in the 1st with direct neighbours

# initialize main lists
qm = []      # list with QM atoms
mm = []      # list with MM atoms
cent    = [] # metal center
shell01 = [] # direct neighbours
shell02 = [] # rings to the direct neighbours

if bool(True):
  printf("Start the actual QM/MM partioning process\n")

################################################################################
# 01 - find metal centers
################################################################################

qm_cnt = 1
for at in mol.atoms:
  if at.is_metal and at.atomic_number >= 19:
    cent.append(at)
cent = my_clean_list(cent)
qm += cent
qm = my_clean_list(qm)
if bool(True):
  printf("01  - Metal center atom(s)\n")
  my_atom_list(cent, "  ")
  fname = "qm_at%02i.xyz" % qm_cnt
  printf("  Write %s\n", fname)
  list_xyz(qm, fname, "Metal center")

################################################################################
# 02a - find the direct neighbours to the metal center
################################################################################

qm_cnt += 1
for mc in cent:
  for at in mc.neighbours:
    Zneighbors.append(at)
    qm.append(at)
shell01 += Zneighbors
shell01 = my_clean_list(shell01)
qm += Zneighbors
qm = my_clean_list(qm)
Zneighbors = my_clean_list(Zneighbors)
if bool(True):
  printf("02a - Direct neighbours to the netal center(s)\n")
  my_atom_list(Zneighbors, "  ")
  fname = "qm_at%02i.xyz" % qm_cnt
  printf("  Write %s\n", fname)
  list_xyz(qm, fname, "Center, direct neighbours")

################################################################################
# 02b - add the atoms of rings containing direct neighbours
################################################################################

qm_cnt += 1
new_at = [] # list for new atoms
for ne in Zneighbors:
  if len(ne.rings) != 0: # Is the atom member of a ring?
    Zcloserings.append(ne.rings)
    for at in ne.rings[0].atoms: # focus on the smallest ring
      new_at.append(at)
shell01 += new_at
shell01  = my_clean_list(shell01)
qm += new_at
qm  = my_clean_list(qm)
if bool(True):
  printf("02b - Rings containing direct neighbour atoms\n")
  printf("  %i small ring(s) found neighbouring the metal center(s)\n", len(Zcloserings))
  if len(Zcloserings) > 0:
    for ri in Zcloserings:
      print(" ", ri[0])
    printf("02b - Newly added atoms\n")
    my_atom_list(new_at, "  ")
  fname = "qm_at%02i.xyz" % qm_cnt
  printf("  Write %s\n", fname)
  list_xyz(qm, fname, "Center, direct neighbour rings")

################################################################################
# 02c - looking for 1 atom links between the rings
################################################################################

qm_cnt += 1
new_at = [] # list for new atoms
SearchSmarts="[R]~*~[R]"
if bool(True):
  printf("02c - %s Bridges in shell01\n", SearchSmarts)
Bridges = []
# Find linking atoms between rings
Bridges = FindSmarts(SearchSmarts, mol)
if bool(True):
  printf("  %i Bridge(s) in molecule\n", len(Bridges))
# clean up list
del_bridges = []
for i, br in enumerate(Bridges):
  # printf("  %i", i)
  rmbf = bool(False)
  if br[0] not in qm:  # remove bridges not starting with a QM atom
    # printf("  H0")
    rmbf=bool(True)
  if br[2] not in qm:  # remove bridges not ending with a QM atom
    # printf("  H2")
    rmbf=bool(True)
  if br[1] in qm:      # remove bridges with the center atom in the QM list
    # printf("  H1")
    rmbf=bool(True)
  # print()
  if rmbf:
    del_bridges.append(i)
# print(Bridges)
for ind in range(len(del_bridges)-1, -1, -1):
  del Bridges[ind]
  # print(Bridges)
if bool(True):
  printf("  %i Bridge(s) linked with shell01\n", len(Bridges))

# add new atoms
for br in Bridges:
  new_at.append(br[1])
# add atoms
shell01 += new_at
shell01  = my_clean_list(shell01)
qm += new_at
qm  = my_clean_list(qm)
if bool(True):
  if len(Bridges) > 0:
    for br in Bridges:
      for at in br:
        MyLabel = "%s-%i" % (at.atomic_symbol, at.index)
        printf("  %-6s", MyLabel)
      printf("  |")
      NewLinkIndex = InMolBonds(mol, br[0].index, br[1].index)
      if NewLinkIndex >= 0:
        printf("  %3i", NewLinkIndex)
        printf("  %6.4f A", mol.bonds[NewLinkIndex].length)
        printf("  %4.2f A", mol.bonds[NewLinkIndex].ideal_bond_length)
      else:
        printf("Link not in the list of accepted bond!\n")
        exit()
      NewLinkIndex = InMolBonds(mol, br[1].index, br[2].index)
      if NewLinkIndex >= 0:
        printf("  %3i", NewLinkIndex)
        printf("  %6.4f A", mol.bonds[NewLinkIndex].length)
        printf("  %4.2f A", mol.bonds[NewLinkIndex].ideal_bond_length)
      else:
        printf("Link not in the list of accepted bond!\n")
        exit()
      printf("\n")
  printf("  %i atom(s) added to shell01\n", len(new_at))
  if len(new_at) > 0:
    my_atom_list(new_at, "  ")
  fname = "qm_at%02i.xyz" % qm_cnt
  printf("  Write %s\n", fname)
  list_xyz(qm, fname, "Center, shell01 complete")

################################################################################
# 02d - looking for 2 atom links between the rings
################################################################################

qm_cnt += 1
new_at = [] # list for new atoms
SearchSmarts="[R]~*~*~[R]"
if bool(True):
  printf("02c - %s Bridges in shell01\n", SearchSmarts)
  
exit()
  
Bridges = []
# Find linking atoms between rings
Bridges = FindSmarts(SearchSmarts, mol)
if bool(True):
  printf("  %i Bridge(s) in molecule\n", len(Bridges))
# clean up list
del_bridges = []
for i, br in enumerate(Bridges):
  # printf("  %i", i)
  rmbf = bool(False)
  if br[0] not in qm:  # remove bridges not starting with a QM atom
    # printf("  H0")
    rmbf=bool(True)
  if br[2] not in qm:  # remove bridges not ending with a QM atom
    # printf("  H2")
    rmbf=bool(True)
  if br[1] in qm:      # remove bridges with the center atom in the QM list
    # printf("  H1")
    rmbf=bool(True)
  # print()
  if rmbf:
    del_bridges.append(i)
# print(Bridges)
for ind in range(len(del_bridges)-1, -1, -1):
  del Bridges[ind]
  # print(Bridges)
if bool(True):
  printf("  %i Bridge(s) linked with shell01\n", len(Bridges))

# add new atoms
for br in Bridges:
  new_at.append(br[1])
# add atoms
shell01 += new_at
shell01  = my_clean_list(shell01)
qm += new_at
qm  = my_clean_list(qm)
if bool(True):
  if len(Bridges) > 0:
    for br in Bridges:
      for at in br:
        MyLabel = "%s-%i" % (at.atomic_symbol, at.index)
        printf("  %-6s", MyLabel)
      printf("  |")
      NewLinkIndex = InMolBonds(mol, br[0].index, br[1].index)
      if NewLinkIndex >= 0:
        printf("  %3i", NewLinkIndex)
        printf("  %6.4f A", mol.bonds[NewLinkIndex].length)
        printf("  %4.2f A", mol.bonds[NewLinkIndex].ideal_bond_length)
      else:
        printf("Link not in the list of accepted bond!\n")
        exit()
      NewLinkIndex = InMolBonds(mol, br[1].index, br[2].index)
      if NewLinkIndex >= 0:
        printf("  %3i", NewLinkIndex)
        printf("  %6.4f A", mol.bonds[NewLinkIndex].length)
        printf("  %4.2f A", mol.bonds[NewLinkIndex].ideal_bond_length)
      else:
        printf("Link not in the list of accepted bond!\n")
        exit()
      printf("\n")
  printf("  %i atom(s) added to shell01\n", len(new_at))
  if len(new_at) > 0:
    my_atom_list(new_at, "  ")
  fname = "qm_at%02i.xyz" % qm_cnt
  printf("  Write %s\n", fname)
  list_xyz(qm, fname, "Center, shell01 complete")

exit()

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