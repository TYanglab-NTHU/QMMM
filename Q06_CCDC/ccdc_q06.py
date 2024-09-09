# 17.08.2024
#   - Testing bimetallic molecule from Ting Yee
#   - Read precisely one molecule and start to analyze it
# 18.08.2034
#   - Add a serach for briges between first layer rings
#   - Start reorganizing the code to preserve intermediate finds
# 19.08.2024
#   - Add aliphatic bridges between atoms in thelayer
# 20.08.2024
#   - Move to ccdc_q06.py
#   - Look for two atom links between rings
#     -> scope problems with the SMARTS string fixed
# 21.08.2024
#   - Clean code and move functions to the top of the code
#   - Organize the shells in one list
#     * shell[0] contains the metal center for compability reasons
#     * shell[1] direct neighbours and their rings
#     * shell[2] atoms from aromatic rings joined to shell[1]
#       this is an iterative process growing from the center out
#   - verbosity level set by command line parameter
#     some minor points need to be fixed later
# 22.08.2024
#   - Add input file management to the command line for easy testing
#   - Create artificial system to test iterarive ring growth
#   - move the search for short links between rings to the end
#   - iterative search for aromatic rings connected to the quantum core
# 24.08.22
#  - add ring size limit

# Load modules for Python coding
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import subprocess               # call the shell

# Load modules and prepare to read the CCS data base
from ccdc          import io                    # read crystal structures
from ccdc.search   import TextNumericSearch     # search by text or numbers
from ccdc.search   import SubstructureSearch    # search for substructures
import ccdc.search

from ccdc.molecule import Molecule, Atom, Bond  # build a molecule
csd_reader = io.EntryReader('CSD')

################################################################################
# Define global variables, try to limit the scope of remaining variables       #
################################################################################

# file names and related
inp_name = "" # name of the input file
qm_cnt   = 0  # counter for xyz files with QM atoms

# initialize main lists
qm    = [] # list with QM atoms
mm    = [] # list with MM atoms
shell = [] # direct neighbours, shells need to be initialized carefully

# lists with intermediate results
# the Z indicates lists with intermediate results
cent        = []  # metal center atoms
Zneighbors  = []  # direct neighbours to the metal centers
Zcloserings = []  # rings in the 1st with direct neighbours
Zaromrings  = []  # aromatic rings for the iterative search

# flags to control the program flow
VerboseFlag = 1 # flag to control the output level
                # 0 no extra output
                # 1 standard text output & xyz files for the steps
                # 2 extra output for debugging

################################################################################
# fake printf and fprintf as used in C                                         #
################################################################################

def printf(format, *args):
  sys.stdout.write(format % args)
def fprintf(stream, format, *args):
  stream.write(format % args)

################################################################################
# writing the xyz file from molecule                                           #
# input  mol    CSD API molecule containing the geometry information           #
#        name   file name (with suffix)                                        #
#        remark 2n line (comment) in the xyz file                              #
# output None                                                                  #
################################################################################

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

################################################################################
# writing the xyz file from a list of atoms                                    #
# input  alist  list with CSD API atoms containing the geometry information    #
#        name   file name (with suffix)                                        #
#        remark 2nd line (comment) in the xyz file                             #
# output None                                                                  #
################################################################################

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

################################################################################
# write a list of atoms to stdout                                              #
# input  ml  list with CSD API atoms to be listed on stdout                    #
#        indent  opt. arg. to define the beginning of the output lines         #
# output None                                                                  #
################################################################################

def my_atom_list(ml, indent=""):
  for at in ml:
    printf("%s%3i  ",  indent, at.index)
    printf("%2s  ",    at.atomic_symbol)
    printf("%10.6f  ", at.coordinates.x)
    printf("%10.6f  ", at.coordinates.y)
    printf("%10.6f\n", at.coordinates.z)
  return()

################################################################################
# these two come together to sort and clean up (unique elements) lists         #
# my_sort_by_index - sort an atom list by the index of the atom                #
#   the function to process the search key in embedded in the function         #
#   input  ml  list to be sorted                                               #
#   output ml  returns the sorted list                                         #
# my_clean_list - making the list elements unique and sort by index            #
#   input  ml  list to be cleaned up                                           #
#   output ml  returns the cleaned up list                                     #
################################################################################

# sorting
def my_sort_by_index(ml):
  def my_key(ma):
    return(ma.index)
  ml.sort(key=my_key)
  return(ml)

# clean up
def my_clean_list(ml):
  ml=list(set(ml))     # remove dublicates
  ml=my_sort_by_index(ml) # sort the atom list by atom index
  return(ml)

################################################################################
# load a xyz file                                                              #
# input  name  full name of the xyz file                                       #
#        ident identfier, used here to hold a commnent about the molecule      #
# output CSD API molecyle generated from the xyz file                          #
#        A xyz file does not carry the necessary Sybyl information needed to   #
#        classify bonds and rings. Reading a xyz file isnoy suitable for the   #
#        QM/MM project. It is better to convert a xyx file into a mol2 file    #
#        for this project.                                                     #
################################################################################

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
  return(mol)

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
  # create a substructure form the smarts string
  substrcuct = ccdc.search.SMARTSSubstructure(MySmarts)
  # search for the substructure
  substructure_search = SubstructureSearch()
  substructure_search.add_substructure(substrcuct)
  # serach for the substructure
  Matches = substructure_search.search(mol)
  
  # evalualte possible matches and build atom list
  if Matches:
    MyAtoms = []
    for MyMatch in Matches:
      AtomMatchList = []
      for at in MyMatch.match_atoms():
        AtomMatchList.append(at)
      MyAtoms.append(AtomMatchList)
  else:
    MyAtoms = []
  # return list with the results
  return(MyAtoms)

################################################################################
# function my_equal_rings to check if two rings are the same                   #
# input  ring01  1st CSD ring ro be compared                                   #
#        ring02  2nd CSD ring ro be compared                                   #
# output True  the atoms joined by the rings are the same                      #
#        False the atoms joined by the rings are different                     #
################################################################################

def my_equal_rings(ring01, ring02):
  set01 = set(ring01.atoms)
  set02 = set(ring02.atoms)
  if set01 == set02:
    return(bool(True))
  else:
    return(bool(False))

################################################################################
################################################################################
##                                                                            ##
## main code of the test                                                      ##
##                                                                            ##
################################################################################
################################################################################

################################################################################
# parse comand line parameter (this should do foe now)                         #
################################################################################

# get command line parameter
my_arglist = sys.argv.copy() # get command line argument
my_arglist.pop(0) # remove the command name

# look for optional verbose parameter
argpos = -1
for argcnt, my_arg in enumerate(my_arglist):
  if my_arg == "-v":
    argpos = argcnt
    break
my_arglist.pop(argpos)
my_arg = my_arglist.pop(argpos)
if my_arg in ["0", "1", "2"]:
  VerboseFlag = int(my_arg)
else:
  printf("bad command line parameter\n")
  printf("use %s [-v 0,1,2] file name\n", sys.argv[0])
  exit()

# the last surviving argument is the file name
if len(my_arglist) != 1:
    printf("bad command line parameter\n")
    printf("use %s [-v 0,1,2] file name\n", sys.argv[0])
    exit()
inp_name = my_arglist[0]

if VerboseFlag > 0:
  printf("Summarize command line parameter\n")
  printf("  VerboseFlag = %i", VerboseFlag)
  if VerboseFlag == 1:
    printf(" (default or cmd line)\n")
  else:
    printf("\n")
  printf("  Input File  = %s\n", inp_name)

del my_arg, my_arglist, argpos, argcnt

################################################################################
# manage old files and the geometry input file                                 #
################################################################################

# clean up old files
if VerboseFlag > 0:
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

# split the command line argument for the input file
root, ext = os.path.splitext(inp_name)
path, file = os.path.split(inp_name)
path += '/'
base = root.replace(path, '')
if VerboseFlag == 2:
  print("ext  ", ext)
  print("path ", path)
  print("base ", base)

# check for file
if not os.path.isfile(inp_name):
  printf("file %s not found\n", inp_name)
  exit()
# handle the input file by extention
if VerboseFlag > 0: printf("Read input file %s\n", inp_name)
input_type = ext.lower()
if input_type==".xyz":
  if VerboseFlag > 0: printf("  xyz input file (depreciated)\n")
  if bool(False): # old code
    mol = load_xyz(inp_name, "xyz file from disc")
    write_xyz(mol, "back.xyz", "write back "+inp_name+" file")
  else: # improved code, mol2 files are much easier to handle
    mol2_inp_name = inp_name.replace(ext, ".mol2")
    if VerboseFlag > 0: printf("  looking for %s\n", mol2_inp_name)
    if not os.path.isfile(mol2_inp_name):
      if VerboseFlag > 0: printf("  create %s\n", mol2_inp_name)
      try:
        subprocess.run(
          ["obabel", inp_name, "-O", mol2_inp_name],
          stdout = subprocess.DEVNULL,  # redirect stdout
          stderr = subprocess.DEVNULL)  # redirect stderr (used by openbabel)
      except:
        printf("  conversion %s -> %s failed\n", inp_name, mol2_inp_name)
        exit()
    if VerboseFlag > 0: printf("  read %s (better file type)\n", mol2_inp_name)
    mol_reader = io.MoleculeReader(mol2_inp_name)
    mol = mol_reader[0]
    if VerboseFlag > 0: printf("  standarize bonds\n")
    mol.assign_bond_types()
    mol.standardise_aromatic_bonds()
    mol.standardise_delocalised_bonds()
elif input_type==".mol2":
  if VerboseFlag > 0:
    printf("  mol2 input file (preferred)\n")
    printf("  read %s directly\n", inp_name)
  mol_reader = io.MoleculeReader(inp_name)
  mol = mol_reader[0]
  if VerboseFlag > 0: printf("  standarize bonds\n")
  mol.assign_bond_types()
  mol.standardise_aromatic_bonds()
  mol.standardise_delocalised_bonds()
else:
  printf("Wrong input file type\n")
  exit()
# clean up variables not to be used later
del root, file, path, base, ext, mol2_inp_name

################################################################################
# a simple test for aromatic rings, exit code after test                       #
################################################################################

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
smarts_list.append("ccccccC#N")   # benzonitrile
smarts_list.append("cCl")         # Cl-C_ar (pos. test)
smarts_list.append("N#N")         # N2 gas  (neg. test)
smarts_list.append("[Fe][N]")     # Fe-N contacts
smarts_list.append("[Fe][O]")     # Fe-N contacts

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
    printf("  Found %i %s group(s) in the molecule.\n", len(matches), smarts)
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
    printf("  No %s group(s) found.\n", smarts)
printf("  %i functional group(s) found.\n", len(func_groups))

################################################################################
# the actual separation code                                                   #
################################################################################

if VerboseFlag > 0:
  printf("Start the actual QM/MM partioning process\n")

################################################################################
# 01 - find metal centers                                                      #
################################################################################

qm_cnt = 1
shell.append([]) # create shell[0]
for at in mol.atoms:
  if at.is_metal and at.atomic_number >= 19:
    cent.append(at)
cent = my_clean_list(cent)
qm += cent
qm = my_clean_list(qm)
shell[0] = cent.copy()
shell[0] = my_clean_list(shell[0])

if VerboseFlag > 0:
  printf("01  - Metal center atom(s)\n")
  my_atom_list(cent, "  ")
  fname = "qm_at%02i.xyz" % qm_cnt
  printf("  Write %s\n", fname)
  list_xyz(qm, fname, "Metal center")

################################################################################
# 02a - find the direct neighbours to the metal center                         #
################################################################################

qm_cnt += 1
shell.append([]) # create shell[1]
for mc in cent:
  for at in mc.neighbours:
    Zneighbors.append(at)
    qm.append(at)
shell[1] += Zneighbors
shell[1] = my_clean_list(shell[1])
qm += Zneighbors
qm = my_clean_list(qm)
Zneighbors = my_clean_list(Zneighbors)
if VerboseFlag > 0:
  printf("02a - Direct neighbours to the netal center(s)\n")
  printf("  %i direct neighbours found\n", len(Zneighbors))
  if len(Zneighbors) > 0:
    my_atom_list(Zneighbors, "  ")
  fname = "qm_at%02i.xyz" % qm_cnt
  printf("  Write %s\n", fname)
  list_xyz(qm, fname, "Center, direct neighbours")

################################################################################
# 02b - add the atoms of rings containing direct neighbours                    #
################################################################################

qm_cnt += 1
new_at = [] # list for new atoms
for ne in Zneighbors:
  if len(ne.rings)>0 and len(ne.rings)<8 : # Is the atom member of a ring?
    Zcloserings.append(ne.rings)
    for at in ne.rings[0].atoms: # focus on the smallest ring
      new_at.append(at)
shell[1] += new_at
shell[1]  = my_clean_list(shell[1])
qm += new_at
qm  = my_clean_list(qm)
if VerboseFlag > 0:
  printf("02b - Rings containing direct neighbour atoms\n")
  printf("  %i small ring(s) found neighbouring the metal center(s)\n", len(Zcloserings))
  if len(Zcloserings) > 0:
    for ri in Zcloserings:
      if VerboseFlag == 2: print(" ", ri[0])
    printf("02b - Newly added atoms\n")
    my_atom_list(new_at, "  ")
  fname = "qm_at%02i.xyz" % qm_cnt
  printf("  Write %s\n", fname)
  list_xyz(qm, fname, "Center, direct neighbour rings")

################################################################################
# 03 - iteratively add aromatic rings sharing atoms with establishe QM core    #
################################################################################

# first round of aromatic ring joined to inner rings
qm_cnt   += 1     # increase QM file counter
new_at   = []     # list for new atoms
shell.append([])  # new shell[2] for aromatic ring atoms
for sa in shell[1]:
  for ri in sa.rings: # test all rings
    if ri.is_aromatic:
      for at in ri.atoms: # add all aromatic ring atoms not in qm
        if at not in qm:
          new_at.append(at)
new_at = my_clean_list(new_at)
for at in new_at:
  for ri in at.rings:
    if ri.is_aromatic and len(ri) <= 10:
      Zaromrings.append(ri)
Zaromrings = list(set(Zaromrings))
shell[2] = new_at.copy()
shell[2] = my_clean_list(shell[2])
qm += new_at
qm = my_clean_list(qm)
if VerboseFlag > 0:
  printf("03a - aromatic rings having shared atoms shell[1]\n")
  printf("  %i atoms found\n", len(new_at))
  if len(new_at) > 0:
    my_atom_list(new_at, "  ")
  fname = "qm_at%02i.xyz" % qm_cnt
  printf("  write %s with the 0th level aromatic rings atoms\n", fname)
  list_xyz(qm, fname, "Center, aromatic rings")

# start growing the number of joined aromatic ring
qm_cnt   += 1     # increase QM file counter
new_at   = []     # list for new atoms

if VerboseFlag > 0:
  printf("03b - iterative search for joined aromatic rings\n")
  printf("  %i rings passed from the previous step\n", len(Zaromrings))
if VerboseFlag == 2:
  for ri in Zaromrings:
    for at in ri.atoms:
      printf("  %3i", at.index)
    print()
if VerboseFlag > 0:
  printf("  cleaning the list of rings\n")
del_ring = []
for i in range(0, len(Zaromrings)-1, 1):
  for j in range(i+1, len(Zaromrings), 1):
    if my_equal_rings(Zaromrings[i], Zaromrings[j]):
      del_ring.append(j)
if VerboseFlag == 2: print(" ", del_ring)
del_ring = list(set(del_ring))
if VerboseFlag == 2: print(" ", del_ring)
for i in range(len(del_ring)-1, -1, -1):
  Zaromrings.pop(del_ring[i])
if VerboseFlag == 2:
  for ri in Zaromrings:
    for at in ri.atoms:
      printf("  %3i", at.index)
    print()
if VerboseFlag > 0:
  printf("  %i rings left after cleaning up\n", len(Zaromrings))

if VerboseFlag > 0:
  printf("  Check the new atoms for aromatic rings\n")
if VerboseFlag == 2:
  for i, ri in enumerate(Zaromrings):
    for at in ri.atoms:
      printf("  %2i", i)
      printf("  %2s", at.atomic_symbol)
      printf("  %3i", at.index)
      printf(" -> ")
      printf(" %2i", len(at.rings))
      for zr in at.rings:
        if zr.is_aromatic:
          printf("  A")
        else:
          printf("  -")
      printf("\n")
for ri in Zaromrings:
  for at in ri.atoms:
    for zr in at.rings:
      if zr.is_aromatic:
        for ra in zr.atoms:
          if ra not in qm:
            new_at.append(ra)
new_at = my_clean_list(new_at)
if VerboseFlag > 0:
  printf("  %i new atoms detected\n", len(new_at))
  if len(new_at) > 0:
    my_atom_list(new_at, "  ")

# clean up variables
del i, j, at, new_at, sa, ri, zr, del_ring, fname

exit()

################################################################################
# 02c - looking for 1 atom links between the rings
################################################################################

qm_cnt += 1
new_at = [] # list for new atoms
SearchSmarts="[R][A][R]"
if VerboseFlag > 0:
  printf("02c - %s Bridges in shell[1]\n", SearchSmarts)
Bridges = []
# Find linking atoms between rings
Bridges = FindSmarts(SearchSmarts, mol)
if VerboseFlag > 0:
  printf("  %i Bridge(s) in molecule\n", len(Bridges))
# clean up list
del_bridges = []
for i, br in enumerate(Bridges):
  if VerboseFlag == 2:  printf("  %i", i)
  rmbf = bool(False)
  if br[0] not in qm:  # remove bridges not starting with a QM atom
    if VerboseFlag == 2: printf("  F")
    rmbf=bool(True)
  else:
    if VerboseFlag == 2: printf("  T")
  if br[1] in qm:      # remove bridges with the center atom in the QM list
    if VerboseFlag == 2: printf("  F")
    rmbf=bool(True)
  else:
    if VerboseFlag == 2: printf("  T")
  if br[2] not in qm:  # remove bridges not ending with a QM atom
    if VerboseFlag == 2: printf("  F")
    rmbf=bool(True)
  else:
    if VerboseFlag == 2: printf("  T")
  if VerboseFlag == 2: print()
  if rmbf:
    del_bridges.append(i)
for ind in range(len(del_bridges)-1, -1, -1):
  del Bridges[ind]
if VerboseFlag > 0:
  printf("  %i Bridge(s) linked with shell[1]\n", len(Bridges))

# add new atoms
for br in Bridges:
  new_at.append(br[1])
# add atoms
shell[1] += new_at
shell[1]  = my_clean_list(shell[1])
qm += new_at
qm  = my_clean_list(qm)
if VerboseFlag > 0:
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
  printf("  %i atom(s) added to shell[1]\n", len(new_at))
  if len(new_at) > 0:
    my_atom_list(new_at, "  ")
  fname = "qm_at%02i.xyz" % qm_cnt
  printf("  Write %s\n", fname)
  list_xyz(qm, fname, "Center, shell[1] complete")
# clean up
del Bridges

################################################################################
# 02d - looking for 2 atom links between the rings
################################################################################

qm_cnt += 1
new_at = [] # list for new atoms
SearchSmarts="[R][A][A][R]"
if VerboseFlag > 0:
  printf("02d - %s Bridges in shell[1]\n", SearchSmarts)
Bridges = []
# Find linking atoms between rings
Bridges = FindSmarts(SearchSmarts, mol)
if VerboseFlag > 0:
  printf("  %i Bridge(s) in molecule\n", len(Bridges))
# clean up list
del_bridges = []
for i, br in enumerate(Bridges):
  if VerboseFlag == 2:
    printf("  %i  ", i)
    print(br, end="")
  rmbf = bool(False)
  if br[0] not in qm:  # remove bridges not starting with a QM atom
    if VerboseFlag == 2: printf("  F")
    rmbf=bool(True)
  else:
    if VerboseFlag == 2: printf("  T")
  if br[1] in qm:      # remove bridges with the center atom in the QM list
    if VerboseFlag == 2: printf("  F")
    rmbf=bool(True)
  else:
    if VerboseFlag == 2: printf("  T")
  if br[2] in qm:      # remove bridges with the center atom in the QM list
    if VerboseFlag == 2: printf("  F")
    rmbf=bool(True)
  else:
    if VerboseFlag == 2: printf("  T")
  if br[3] not in qm:  # remove bridges not ending with a QM atom
    if VerboseFlag == 2: printf("  F")
    rmbf=bool(True)
  else:
    if VerboseFlag == 2: printf("  T")
  if VerboseFlag == 2: print()
  if rmbf:
    del_bridges.append(i)
# print(Bridges)
for ind in range(len(del_bridges)-1, -1, -1):
  del Bridges[ind]
  # print(Bridges)
if VerboseFlag > 0:
  printf("  %i Bridge(s) linked with shell[1]\n", len(Bridges))
# add new atoms
for br in Bridges:
  new_at.append(br[1])
# add atoms
shell[1] += new_at
shell[1]  = my_clean_list(shell[1])
qm += new_at
qm  = my_clean_list(qm)
if VerboseFlag > 0:
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
  printf("  %i atom(s) added to shell[1]\n", len(new_at))
  if len(new_at) > 0:
    my_atom_list(new_at, "  ")
  fname = "qm_at%02i.xyz" % qm_cnt
  printf("  Write %s\n", fname)
  list_xyz(qm, fname, "Center, shell[1] complete")
# clean up
del Bridges

exit()

# the search for ligand rings ends here
# the search for ligand chains begins
# add functional groups

qm_cnt += 1
shell03=[]
for ne_at in shell[1]: # Loop over direct neighbors
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
