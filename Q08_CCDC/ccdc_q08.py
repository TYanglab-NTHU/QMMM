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
#   - begin iterative search for aromatic rings connected to the quantum core
# 24.08.2024
#  - add ring size limit
#  - move to ccdc_q07.py and keep ccdc_q06.py as backup
# 25.08.2024
#  - the order of atoms in the SMARTS string is not reflected in the result
#    -> rewrite the code for 1 atom bridges
# 26.08.2024
#  - the order of atoms in the SMARTS string is not reflected in the result
#    -> rewrite the code for 2 atom bridges
#    -> add distace parameter for the neighbours
# 27.08.2024
#  - move to ccdc_q08.py
#  - start with the compartmentalization of the Python code
#    * limit the scope of variables :)
#    * make it easiear to import my code into other projects
#    * increase the redability of the main code
#  - Each step gets its on shell to make the code more clear
#    (One atom can show up shells, but only once in the qm and mm lists)
#    -> rewrite shell definitions form 21.08.2024
#       shell[0]  metal center
#       shell[1]  direct neighbours
# 28.08.2024
#  - Continue with the clean up
#    -> I rewrote the iterative search for aromatics rings
#       The code is cleaner, but appears to be slower (?)
#    -> New shells
#       shell[2]  distant neighbours, 80% VdW radii sum
#       shell[3]  rings containing direct neighbours
#       shell[4]  1st round of aromatic rings
#       shell[5]  outer aromatic rings, iterative search
#       shell[6]  1 atom links between QM atoms
#       shell[7]  2 atom links between QM atoms
#       shell[8]  H atoms directly attached to the QM atoms
#       shell[9]  direct neighbours
#  - Nove to ccdc_q09.py


# Load modules for Python coding
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import subprocess               # call the shell
import math                     # mathematical  functions

# Load modules and prepare to read the CCS data base
from ccdc          import io                    # read crystal structures
from ccdc.search   import TextNumericSearch     # search by text or numbers
from ccdc.search   import SubstructureSearch    # search for substructures
import ccdc.search

from ccdc.molecule import Molecule              # Build a molecule
from ccdc.molecule import Atom                  # Atomic data and properties
from ccdc.molecule import Bond                  # Bond properties
csd_reader = io.EntryReader('CSD')

################################################################################
# Define global variables and try to limit the scope of remaining variables    #
################################################################################

# file names and related
inp_name = "" # name of the input file
qm_cnt   = 0  # counter for xyz files with QM atoms
stp_cnt  = 0  # general step counter

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
# general purpose function                                                     #
# fake printf and fprintf as used in C                                         #
################################################################################

def printf(format, *args):
  sys.stdout.write(format % args)
def fprintf(stream, format, *args):
  stream.write(format % args)

################################################################################
# general purpose function                                                     #
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
# general purpose function                                                     #
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
# general purpose function                                                     #
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
# general purpose function                                                     #
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
# functions to make it easier with the API                                     #
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
# functions to make it easier with the API                                     #
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
# functions to make it easier with the API                                     #
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
# functions to make it easier with the API                                     #
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
# functions to make it easier with the API                                     #
# calculate distance between atoms, no bond needed                             #
# input  at1  1st CSD atom                                                     #
#        at2  2nd CSD atom                                                     #
# output dist float with the distance between the atoms                        #
#             (The unit depends on the unit of the cartesian coordinates)      #
################################################################################

def atom_dist(at1, at2):
  dist = 0.0
  dist += (a1.coordinates.x-a2.coordinates.x)**2
  dist += (a1.coordinates.y-a2.coordinates.y)**2
  dist += (a1.coordinates.z-a2.coordinates.z)**2
  dist = math.sqrt(dist)
  return(dist)

################################################################################
# functions to for the compartimensatiom of the  QM/MM separation              #
# this function summarizes an ibdividual step of the search                    #
#   global var VerboseFlag will be used to controll the output                 #
# input  cnt step counter                                                      #
#        nlist new atom list                                                   #
#        flist atom list send to file                                          #
#        mytxt string for the description                                      #
################################################################################

def summarize_step(cnt, nlist, flist, mytxt):
  # is VerboseFlag globally defined and in range?
  if "VerboseFlag" not in globals():
    print("Variable 'VerboseFlag' not globally defined")
    exit()
  if VerboseFlag<0 or VerboseFlag>2:
    print("Control flag 'VerboseFlag' out of range.")
    printf("Expcted 0, 1, or 2 but found %i\n", VerboseFlag)
    exit()
  # output to stdout
  if VerboseFlag == 0:
    return bool(True)
  printf("Step %02i - %s\n", cnt, mytxt)
  printf("  %i new atoms found\n", len(nlist))
  if VerboseFlag == 2 and len(nlist)>0:
    my_atom_list(nlist, "  ")
  # write atom list to disk
  fname = "stp%02i.xyz" % cnt
  printf("  Write %s\n", fname)
  list_xyz(flist, fname, mytxt)
  return bool(True)

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
if VerboseFlag > 0: printf("Clean up old files\n")
for file in glob.glob("qm_at*xyz"):
  os.remove(file)
for file in glob.glob("stp*xyz"):
  os.remove(file)
if os.path.exists("./mm_at00.xyz"): os.remove("./mm_at00.xyz")
if os.path.exists("./combi.xyz"): os.remove("./combi.xyz")

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

stp_cnt = 1
shell.append([]) # create shell[0] for central metal atoms

for at in mol.atoms:
  if at.is_metal and at.atomic_number >= 19:
    shell[0].append(at)
shell[0] = my_clean_list(shell[0])
qm += shell[0]
qm = my_clean_list(qm)
summarize_step(stp_cnt, shell[0], qm, "metal center")

################################################################################
# 02 - find the direct neighbours to the metal center                          #
################################################################################

stp_cnt += 1
shell.append([]) # create shell[1] for direct neighbours

for mc in shell[0]:
  for at in mc.neighbours:
    shell[1].append(at)
    qm.append(at)
shell[1] = my_clean_list(shell[1])
qm = my_clean_list(qm)
summarize_step(stp_cnt, shell[1], qm, "direct neighbours")

################################################################################
# 03 - additional search based on distance, use the VdW radii for validity     #
#      The atom radius seems NOT to be defined in the CSD API. A valid         #
#      alternative is the Van der Waals Radius. If the VdW radius is not       #
#      defined for the CSD API, the API returns a default value of 2.0 A.      #
#      I set the cut-off value arbitrarily at 80% of the VdW sum.              #
#      The pure distance cut-off is set at 3.0 A.                              #
################################################################################

stp_cnt += 1
shell.append([]) # create shell[2] for distance neighbours
bnd_cnt = 0      # counter for new bonds
vdw_sum = 0.0    # sum of the VdW radii of a long distance contact
vdw_fac = 0.8    # cut-off factor for the sum of VdW radii

if VerboseFlag == 2:
  printf("Step 03 - Distant neighbour search\n")
  printf("  <= 3.0 A and <= %.0f%% of the VdW radius\n", vdw_fac*100)
for a1 in shell[0]: # loop over all atoms
  for a2 in mol.atoms: # loop over all atoms
    if a1 == a2 : continue
    dist = atom_dist(a1, a2)
    vdw_sum =  (1.0 * a1.vdw_radius)  # The type of the vdw_property appears to be
    vdw_sum += (1.0 * a2.vdw_radius)  # context dependent. I force it to a float.
    vdw_sum *= vdw_fac                # Apply cut-off factor
    if dist < 3.0 and a2 not in qm:
      if VerboseFlag == 2:
        printf("  %s%i",    a1.atomic_symbol, a1.index)
        printf( "-%s%i",    a2.atomic_symbol, a2.index)
        printf("  %6.4f",   dist)
        printf("  %6.4f",   vdw_sum)
        if dist<=vdw_sum: printf("  Passed\n")
        else: printf("  too long\n")
      if dist<=vdw_sum:
        shell[2].append(a2)
        if InMolBonds(mol, a1, a2) == -1:
          mol.add_bond(1, a1, a2)
          bnd_cnt += 1
shell[2] = my_clean_list(shell[2])
qm += shell[2]
qm = my_clean_list(qm)
summarize_step(stp_cnt, shell[2], qm, "distance based neighbours")

################################################################################
# 04 - add the atoms of rings containing direct neighbours                     #
################################################################################

stp_cnt += 1
shell.append([]) # create shell[3] for distance neighbours

for ne in my_clean_list(shell[1]+shell[2]):
  if len(ne.rings)>0 and len(ne.rings[0])<=10 : # Is the atom member of a ring?
    for at in ne.rings[0].atoms: # focus on the smallest ring
      if at not in qm:
        shell[3].append(at)
shell[3] = my_clean_list(shell[3])
qm += shell[3]
qm  = my_clean_list(qm)
summarize_step(stp_cnt, shell[3], qm, "rings containing neighbours atoms")

################################################################################
# 05 - iteratively add aromatic rings sharing atoms with establishe QM core    #
#      first round of aromatic ring joined to inner rings                      #
################################################################################

stp_cnt   += 1    # increase QM file counter
shell.append([])  # create shell[4] for inner aromatic rings

# create an intermediate test list for atoms to be tested
itl = my_clean_list(shell[1]+shell[2]+shell[3])
# find new MM aromatic ring atoms
for sa in itl:
  for ri in sa.rings: # test all rings
    if ri.is_aromatic:
      for at in ri.atoms: # add all aromatic ring atoms not in qm
        if at not in qm:
          shell[4].append(at)
shell[4] = my_clean_list(shell[4])
qm += shell[4]
qm = my_clean_list(qm)
summarize_step(stp_cnt, shell[4], qm, "central aromatic rings")

# clean up variables
del itl

################################################################################
# 06 - grow the number of joined aromatic ring                                 #
################################################################################

stp_cnt   += 1      # increase QM file counter
shell.append([])    # create shell[5] for outer aromatic rings
it_cnt = 0          # counter for iterative rounds          

# create an intermediate test list for atoms to be tested
itl = my_clean_list(shell[1]+shell[2]+shell[3]+shell[4])

while bool(True):   # infinite loop for the iterative search
  new_at = []       # atoms found by each loop of the iterative search
  it_cnt += 1       # increase counter
  # print current step to stdout
  if VerboseFlag == 2:
    printf("Step 06 - iterative search for aromatic rings\n")
    printf("  %i atoms found in the previous step\n", len(new_at))
  # check whether the new atoms are part of an aromatic ring and that these atoms
  # are not in the QM core
  for at in itl:
    for ri in at.rings:
      if ri.is_aromatic:
        for tat in ri.atoms:
          if tat not in qm:
            new_at.append(tat)
  if len(new_at) == 0:
    break
  else:
    shell[5] += new_at
    shell[5] = my_clean_list(shell[5])
    itl = my_clean_list(shell[1]+shell[2]+shell[3]+shell[4]+shell[5])
    new_at = []
    qm += shell[5]
    qm = my_clean_list(qm)
if VerboseFlag == 2:
  printf("  %i iterative search rounds needed\n", it_cnt)
summarize_step(stp_cnt, shell[4], qm, "iterative search for aromatic rings")

# clean up variables
del itl, it_cnt

################################################################################
# 07 - looking for 1 atom links between QM atoms                               #
################################################################################

stp_cnt   += 1      # increase QM file counter
shell.append([])    # create shell[6] for 1 atom links

if VerboseFlag == 2:
  printf("Step %02i - 1 atom links", stp_cnt-1)
  printf("  list of all found links\n")
for sa in qm:
  for ea in qm:
    if sa.index == ea.index:
      continue
    for la in sa.neighbours:
      if la in ea.neighbours and la not in qm:
        shell[6].append(la)
        if VerboseFlag == 2:
          printf("  %s%i-", sa.atomic_symbol, sa.index)
          printf("%s%i-", la.atomic_symbol, la.index)
          printf("%s%i\n", ea.atomic_symbol, ea.index)
shell[6] = my_clean_list(shell[6])
qm += shell[6]
qm = my_clean_list(qm)
summarize_step(stp_cnt, shell[6], qm, "1 atom links")

################################################################################
# 08 - looking for 2 atom links between the rings                              #
################################################################################

stp_cnt   += 1      # increase QM file counter
shell.append([])    # create shell[7] for 2 atom links

if VerboseFlag == 2:
  printf("Step %02i - 2 atom links", stp_cnt-1)
  printf("  list of all found links\n")

for b0 in qm:
  for b1 in b0.neighbours:
    for b2 in b1.neighbours:
      for b3 in b2.neighbours:
        if VerboseFlag == 2:
          printf("  %s%i-", b0.atomic_symbol, b0.index)
          printf("%s%i-",   b1.atomic_symbol, b1.index)
          printf("%s%i-",   b2.atomic_symbol, b2.index)
          printf("%s%i  ",  b3.atomic_symbol, b3.index)
          if (b1 not in qm and b2 not in qm and b3 in qm):
            printf("Add\n")
          else:
            printf("Del\n")
        if (b1 not in qm and b2 not in qm and b3 in qm):
          shell[7].append(b1)
          shell[7].append(b2)
shell[7] = my_clean_list(new_at)
qm += shell[7]
qm = my_clean_list(qm)
summarize_step(stp_cnt, shell[7], qm, "2 atom links")

################################################################################
# 09 - adding single hydrogen atoms                                            #
################################################################################

stp_cnt   += 1      # increase QM file counter
shell.append([])    # create shell[8] for H atoms links

for at in qm:
  for na in at.neighbours:
    if na.atomic_number == 1 and na not in qm:
      shell[8].append(na)
shell[8] = my_clean_list(shell[8])
qm += shell[8]
qm = my_clean_list(qm)
summarize_step(stp_cnt, shell[8], qm, "H atoms")

################################################################################
# Final step - create MM layer
#   any atom not QM is MM
#   combine both sets in one xyz file for output
################################################################################

stp_cnt   += 1      # increase file counter
mm = []             # create the list for MM atoms

# create mm layer
for at in mol.atoms:
  if at in qm:
    pass
  else:
    mm.append(at)
mm = my_clean_list(mm)
summarize_step(stp_cnt, mm, mm, "only MM atoms")

# create combined output file
combi=[]
if VerboseFlag > 0:
  printf("Combined output:\n")
  printf("  QM atoms: %i\n", len(qm))
  printf("  MM atoms: %i\n", len(mm))
combi = qm + mm
fname = "combi.xyz"
comment = "QM: %i atoms, MM: %i atoms" % (len(qm), len(mm))
list_xyz(combi, fname, comment)

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
