# 17.08.2024
# - Testing bimetallic molecule from Ting Yee
# - Read precisely one molecule and start to analyze it
# 18.08.2034
# - Add a serach for briges between first layer rings
# - Start reorganizing the code to preserve intermediate finds
# 19.08.2024
# - Add aliphatic bridges between atoms in thelayer
# 20.08.2024
# - Move to ccdc_q06.py
# - Look for two atom links between rings
#   -> scope problems with the SMARTS string fixed
# 21.08.2024
# - Clean code and move functions to the top of the code
# - Organize the shells in one list
#   * shell[0] contains the metal center for compability reasons
#   * shell[1] direct neighbours and their rings
#   * shell[2] atoms from aromatic rings joined to shell[1]
#     this is an iterative process growing from the center out
# - verbosity level set by command line parameter
#   some minor points need to be fixed later
# 22.08.2024
# - Add input file management to the command line for easy testing
# - Create artificial system to test iterarive ring growth
# - move the search for short links between rings to the end
# - begin iterative search for aromatic rings connected to the quantum core
# 24.08.2024
# - add ring size limit
# - move to ccdc_q07.py and keep ccdc_q06.py as backup
# 25.08.2024
# - the order of atoms in the SMARTS string is not reflected in the result
#   -> rewrite the code for 1 atom bridges
# 26.08.2024
# - the order of atoms in the SMARTS string is not reflected in the result
#   -> rewrite the code for 2 atom bridges
#   -> add distace parameter for the neighbours
# 27.08.2024
# - move to ccdc_q08.py
# - start with the compartmentalization of the Python code
#   * limit the scope of variables :)
#   * make it easiear to import my code into other projects
#   * increase the redability of the main code
# - Each step gets its on shell to make the code more clear
#   (One atom can show up shells, but only once in the qm and mm lists)
#   -> rewrite shell definitions form 21.08.2024
#      shell[0]  metal center
#      shell[1]  direct neighbours
# 28.08.2024
# - Continue with the clean up
#   -> I rewrote the iterative search for aromatics rings
#      The code is cleaner, but appears to be slower (?)
#   -> New shells
#      shell[2]  distant neighbours, 80% VdW radii sum
#      shell[3]  rings containing direct neighbours
#      shell[4]  1st round of aromatic rings
#      shell[5]  outer aromatic rings, iterative search
#      shell[6]  1 atom links between QM atoms
#      shell[7]  2 atom links between QM atoms
#      shell[8]  H atoms directly attached to the QM atoms
#      shell[9]  all MM atoms, no QM atoms
# - Nove to ccdc_q09.py
# - remove the quick test for aromatic and fully conjugated rings
#   Return to ccdc_q08.py, if I needed
# - started with con_unit
#   Can the 2 atoms in question part of mesomeric chain?
# 04.09.2024
# - Remove the addition of new bonds for distant neighbours
#   to preserve the CSD assignment of small rings.
# 05.09.2024
# - move to ccdc_q10.py
# - general function added to remove a list of global variables
#   to limit their scope and clean up the code
# - I do thisto build a paralelle OpenBabel object holding a
#   copy of the CSD molecule. If something goes wrong, go back
#   to ccdc_q09 and write additional code to get the chemical
#   properties of an atom.
#   * The OpenBabel molecule object will part of preocessing the 
#     input file.
#   * Python module OpenBabel used to convert a xyzfile to mol2
#     and write the resulting mol2 file to disk. The OpenBabel
#     objects are deleted after the conversion
#   * ObenBabel object implemented. I try to ensure consistency
#     by reading the same mol2 file both tests. An additional
#     test makes sure that both geometries are exactly the same.
# 06.09.24
# - Use OpenBabel in the detection of conjugated triads
# - Grow conjugated chains
# 07.09.24
# - Finish conjugated chain growth
# 10.09.24
# - Growth of conjugated chains finished
#   Next: Marry the chain code with the detection of rings
# - Move to ccdc_q11.py
#   To marry the search for conjugated chains with th iterative
#   search for aromatic rings, I have to do major changes to
#   the base structure of the code.
#   ATTENTION: The directory entry has to be last action of a step.
#     stp 01 -> cent,    shell[0]  metal center
#     stp 02 -> close,   shell[1]  close neighbors
#     stp 03 -> dist,    shell[2]  distant neighbors
#     stp 04 -> nrings,  shell[3]  rings with neighboring atoms
#     stp 05 -> cArings, shell[4]  central aromatic rings
#     stp 06 -> gArings, shell[5]  central aromatic rings
# - Screw up in the output from step 06 fixed. The search results
#   were oK and did needed to be changed.
# 11.09.24
# - Continue with the work on the dictionary
#     stp 07 -> aLink,   shell[6]  1 atom bridges
#     stp 08 -> aaLink,  shell[7]  2 atom bridges
#     stp 09 -> hatom,   shell[8]  H atoms attached to the QM core
# - Transformation of shells[0] to shell[8] finished and double
#   checked using <grep -n "shell\[8\]" ccdc_q11.py>
# - Make the serach for conjugated units movable
#   The conjugated atoms were collected in the old shell[9]
#     stp xx -> conju,   shell[9]  atoms in conjugated chains
# - Move to ccdc_q12.py
#   To make the individual search steps movable, I want to turn
#   them into functions. This also helps with the scope of the
#   variables. The functions will return a list with the new atoms.
#   The functions will probe the necessary global variables so that
#   they don't need to be passed to the function. [But it might be
#   the cleaner way to go.]
#   * conjugated chains done
#   * lone H atoms done
#   * 2-atom links
#   * 1-atom links
#   * grow of aromatic rings done
#   * inner aromatic rings done
#   * neighboring rings done
#   * distant neighbors done
#   * close neighbors done
#   * metal center done
# - It looks like moving the search blocks into individual functions
#   removed the need for the shell list. I will delete it from the
#   code and use the dictionary dnts instead.
#   The code doesn't crash after I removed the shell list.
# 12.09.2024
# - Before I merge both search metods, I slim the main code
# - Merge both search strategies
#   Tests with chain01 to chain06 worked out well
# - GitHub commit 240912 - 13:50
# - Cleanup code
#   * command line parser stays
#   * clean up of old files stays
#   * reading input files goes into a function

# Load modules for Python coding
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
# import subprocess               # call the shell
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

# Load molecules to use OpenBabel within the QM/MM separation code
from openbabel import openbabel as ob           # basic OpenBabel

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
dnts  = {} # Dictionary to look up shells by name (Name To Shell)

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
# function to delte global variables to limit the scope of these variables     #
# input  vlist  a list with the names (strings) of the variables to be         #
#               deleted.                                                       #
# output anz    integer with the number of deleted variables                   #
################################################################################

def clean_glob_var(vlist):
  if "VerboseFlag" not in globals():
    print("Variable 'VerboseFlag' not globally defined")
    exit()
  anz = 0
  g = globals()
  for var in vlist:
    if var in globals():
      if VerboseFlag==2: printf("%s in globals\n", var)
      anz += 1
      del g[var]
    else:
      if VerboseFlag==2: printf("%s NOT in globals\n", var)
  if VerboseFlag==2: printf("%i global variables deleted\n", anz)
  return(anz)

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
# general purpose function                                                     #
# create a custom label for an atom                                            #
# input  at     atom as defined by the CSD Api                                 #
# output mlabel custom label for the atom                                      #
################################################################################

def MyLabel(at):
  labelstrg = "%s-%i" % (at.atomic_symbol, at.index)
  return(labelstrg)

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
  dist += (at1.coordinates.x-at2.coordinates.x)**2
  dist += (at1.coordinates.y-at2.coordinates.y)**2
  dist += (at1.coordinates.z-at2.coordinates.z)**2
  dist = math.sqrt(dist)
  return(dist)

################################################################################
# functions to make it easier with the API                                     #
# create a string to label a triad of atoms                                    #
# input  trip a list of 3 atoms (no sanity control)                            #
# output mstr string with the label for the triad                              #
################################################################################

def triad_string(trip):
  mstr  = "%s%i" % (trip[0].atomic_symbol, trip[0].index)
  mstr += "-"
  mstr += "%s%i" % (trip[1].atomic_symbol, trip[1].index)
  mstr += "-"
  mstr += "%s%i" % (trip[2].atomic_symbol, trip[2].index)
  return(mstr)

################################################################################
# functions to to help with the QM/MM separation                               #
# can theese tree atoms be a unit of conjugated chain?                         #
# this code is looking only for standard 2nd row bond types                    #
# A doublebond folled by a single bond is not enough to detect a unit, because #
# all atoms need to  have p-orbitals ready to engage in pi-bonding.            #
# input  mymol molecule containing the atoms                                   #
#        at1   1st CSD atom                                                    #
#        at2   2nd CSD atom                                                    #
#        at3   3rd CSD atom                                                    #
# output utype interger to indicate the type of 3-atom unit                    #
#        -3  triple-single                                                     #
#        -2  double-single                                                     #
#         0  not a viable candidate (default)                                  #
#         2  single-double                                                     #
#         3  single-triple                                                     #
# bond type definitions form                                                   #
# downloads.ccdc.camd.ac.uk/documentation/API/descriptive_docs/molecule.html   #
################################################################################

def con_unit(mymol, at1, at2, at3):
  # No H atoms allowed
  if 1 in [at1.atomic_number, at2.atomic_number, at3.atomic_number]:
    # print("No H atoms allowd")
    return(0)
  if VerboseFlag==2:
    printf("Testing triad %s\n", triad_string([at1, at2, at3]))
  if bool(True):
    # Use the number of neighbours to estimate the hybridization state
    # Each atom shukd have less then 4 direct neighbours
    if len(at1.neighbours)>3 or len(at1.neighbours)<1: return(0)  # not a candidate
    if len(at2.neighbours)>3 or len(at2.neighbours)<1: return(0)  # not a candidate
    if len(at3.neighbours)>3 or len(at3.neighbours)<1: return(0)  # not a candidate
  else:
    # I try to read the OpenBabel hybridization information 
    # 1 for sp, 2 for sp2, 3 for sp3, 4 for sq. planar, 5 for trig. bipy, 6 for octahedral
    # https://openbabel.org/api/3.0/classOpenBabel_1_1OBAtom.shtml
    printf("  Hyb")
    obatom=obmol.GetAtom(at1.index+1)
    printf("  %i", obatom.GetHyb())
    if obatom.GetHyb() not in [1, 2]:
      return(0)
    obatom=obmol.GetAtom(at2.index+1)
    printf("  %i", obatom.GetHyb())
    if obatom.GetHyb() not in [1, 2]:
      return(0)
    obatom=obmol.GetAtom(at3.index+1)
    printf("  %i\n", obatom.GetHyb())
    if obatom.GetHyb() not in [1, 2]:
      return(0)
  # conjugated chains can start or end at an aromatic rind, but should not be
  # part of an aromatic ring
  obatom=obmol.GetAtom(at2.index+1)
  if obatom.IsAromatic(): return(0)  # not a candidate
  # get the indices of the bonds linking the atoms and check whether the atoms
  # are bonded
  utype = 0 # set default, not a candidate
  ind1 = InMolBonds(mymol, at1.index, at2.index)
  ind2 = InMolBonds(mymol, at2.index, at3.index)
  if VerboseFlag==2:
    printf("  ind1: %i, ind2 %i\n", ind1, ind2)
  if ind1 == -1: return(0) # no bond between at1 and at2
  if ind2 == -1: return(0) # no bond between at2 and at3
  # Reply by bond type
  # How likely isit that I catch bond types 7 and 9 ?
  if VerboseFlag==2:
    printf("  %s", MyLabel(at1))
    printf(" %s", mymol.bonds[ind1].bond_type)
    printf(" %s  ", MyLabel(at2))
    printf(" %s", mymol.bonds[ind2].bond_type)
    printf(" %s\n", MyLabel(at3))
  if mymol.bonds[ind1].bond_type == 3 and mymol.bonds[ind2].bond_type == 1:
    return(-3)
  if mymol.bonds[ind1].bond_type == 2 and mymol.bonds[ind2].bond_type == 1:
    return(-2)
  if mymol.bonds[ind1].bond_type == 1 and mymol.bonds[ind2].bond_type == 3:
    return(3)
  if mymol.bonds[ind1].bond_type == 1 and mymol.bonds[ind2].bond_type == 2:
    return(2)
  return(utype)

################################################################################
# functions to to help with the QM/MM separation                               #
# test wether a conjugated 3-atom unit starts at the given atom                #
# The search does not detct all possible units. The serach stops after the one #
# valid unit has been found.                                                   #
# input  at1        1st atom (CDC atoms type) of the triad                     #
# output con_level  integer descibing the bond pattern of the unit             #
#        triad      list with the 3 atoms of the conjugated unit               #
################################################################################

def check_conjugation(at1):
  triad = []
  for at2 in at1.neighbours:
    #if at2 in qm:
    #  continue
    for at3 in at2.neighbours:
      if at3 in qm:
        continue
      if at1 == at3:
        continue
      con_level = con_unit(mol, at1, at2, at3)
      if con_level != 0:
        triad += [at1, at2, at3]
        return(con_level, triad)
  return(0, triad)

################################################################################
# functions to to help with the QM/MM separation                               #
# find chains of conjugated 3-atom units                                       #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
################################################################################

def conjugated_chains():
  new_atoms = []
  # are the necessary global variables defined?
  gvar_list = ["qm", "VerboseFlag"]
  for gvar in gvar_list:
    if gvar not in globals():
      printf("Variable %s NOT in globals.\n", gvar)
      exit()
    elif VerboseFlag==2: printf("Variable %s in globals.\n", gvar)
  global qm
  # actual search
  while bool(True): # infinite loop to detect chains of units
    cl1 = 0
    trip = []
    gotcha = 0
    for at1 in my_clean_list(qm): # Loop over all QM atoms
      if VerboseFlag==2:
        printf("Testing %s-%i\n", at1.atomic_symbol, at1.index)
      cl1, trip = check_conjugation(at1)
      if cl1 != 0:
        if VerboseFlag==2:
          printf("  %-10s  %2i\n", triad_string(trip), cl1)
        for nat in trip:
          if nat not in qm:
            gotcha = 1
            new_atoms.append(nat)
    new_atoms = my_clean_list(new_atoms)
    qm += new_atoms
    qm = my_clean_list(qm)
    if gotcha == 0:
      break
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# find H atoms attached to the QM core                                         #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
################################################################################

def lone_H_atoms():
  new_atoms = []
  # are the necessary global variables defined?
  gvar_list = ["qm", "VerboseFlag"]
  for gvar in gvar_list:
    if gvar not in globals():
      printf("Variable %s NOT in globals.\n", gvar)
      exit()
    elif VerboseFlag==2: printf("Variable %s in globals.\n", gvar)
  global qm
  # the actual search
  for at in qm:
    for na in at.neighbours:
      if na.atomic_number == 1 and na not in qm:
        new_atoms.append(na)
  new_atoms = my_clean_list(new_atoms)
  qm += new_atoms
  qm = my_clean_list(qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# find two atom bridges                                                        #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
################################################################################

def two_atom_links():
  new_atoms = []
  # are the necessary global variables defined?
  gvar_list = ["qm", "VerboseFlag"]
  for gvar in gvar_list:
    if gvar not in globals():
      printf("Variable %s NOT in globals.\n", gvar)
      exit()
    elif VerboseFlag==2: printf("Variable %s in globals.\n", gvar)
  global qm
  # the actual search
  if VerboseFlag == 2:
    printf("List of all found 2-atom links\n")
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
            new_atoms.append(b1)
            new_atoms.append(b2)
  new_atoms = my_clean_list(new_atoms)
  qm += new_atoms
  qm = my_clean_list(qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# find one atom bridges                                                        #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
################################################################################

def one_atom_links():
  new_atoms = []
  # are the necessary global variables defined?
  gvar_list = ["qm", "VerboseFlag"]
  for gvar in gvar_list:
    if gvar not in globals():
      printf("Variable %s NOT in globals.\n", gvar)
      exit()
    elif VerboseFlag==2: printf("Variable %s in globals.\n", gvar)
  global qm
  # the actual search
  if VerboseFlag == 2:
    printf("List of all found 1-atom links\n")
  for sa in qm:
    for ea in qm:
      if sa.index == ea.index:
        continue
      for la in sa.neighbours:
        if la in ea.neighbours and la not in qm:
          new_atoms.append(la)
          if VerboseFlag == 2:
            printf("  %s%i-", sa.atomic_symbol, sa.index)
            printf("%s%i-", la.atomic_symbol, la.index)
            printf("%s%i\n", ea.atomic_symbol, ea.index)
  new_atoms = my_clean_list(new_atoms)
  qm += new_atoms
  qm = my_clean_list(qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# outward growth of aromatic rings                                             #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def grow_aromatic_rings():
  # are the necessary global variables defined?
  gvar_list = ["qm", "dnts", "VerboseFlag"]
  for gvar in gvar_list:
    if gvar not in globals():
      printf("Variable %s NOT in globals.\n", gvar)
      exit()
    elif VerboseFlag==2: printf("Variable %s in globals.\n", gvar)
  global qm, dnts
  # the actual search
  new_atoms = []  # empty list for serach results
  it_cnt = 0      # counter for iterative rounds          
  # create an intermediate test list for atoms to be tested
  # itl = dnts['close']+dnts['dist']+dnts['nrings']+dnts['cArings']+dnts['conju']
  itl = list(set(qm) - set(dnts['cent']))
  itl = my_clean_list(itl)
  if VerboseFlag == 2:
    printf("Details for the iterative search for aromatic rings\n")
  while bool(True):   # infinite loop for the iterative search
    new_at = []       # atoms found by each loop of the iterative search
    it_cnt += 1       # increase counter
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
      if VerboseFlag == 2:
        printf("  %3i atoms found in step %2i\n", len(new_at), it_cnt)
      new_atoms += new_at
      new_atoms = my_clean_list(new_atoms)
      itl = my_clean_list(itl + new_atoms)
      new_at = []
      qm += new_atoms
      qm = my_clean_list(qm)
  if VerboseFlag == 2:
    printf("  %i iterative search rounds needed\n", it_cnt-1)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# detect the inner aromatic rings                                              #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def inner_aromatic_rings():
  # are the necessary global variables defined?
  gvar_list = ["qm", "dnts", "VerboseFlag"]
  for gvar in gvar_list:
    if gvar not in globals():
      printf("Variable %s NOT in globals.\n", gvar)
      exit()
    elif VerboseFlag==2: printf("Variable %s in globals.\n", gvar)
  global qm, dnts
  # the actual search
  new_atoms = []  # empty list for serach results
  # create an intermediate test list for atoms to be tested
  itl = dnts['close'] + dnts['dist'] + dnts['nrings']
  itl = my_clean_list(itl)
  # find new MM aromatic ring atoms
  for sa in itl:
    for ri in sa.rings: # test all rings
      if ri.is_aromatic:
        for at in ri.atoms: # add all aromatic ring atoms not in qm
          if at not in qm:
            new_atoms.append(at)
  new_atoms = my_clean_list(new_atoms)
  qm += new_atoms
  qm = my_clean_list(qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# find the rings containing neighbors to the metal center                      #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def neighbor_rings():
  # are the necessary global variables defined?
  gvar_list = ["qm", "dnts", "VerboseFlag"]
  for gvar in gvar_list:
    if gvar not in globals():
      printf("Variable %s NOT in globals.\n", gvar)
      exit()
    elif VerboseFlag==2: printf("Variable %s in globals.\n", gvar)
  global qm, dnts
  # the actual search
  new_atoms = []  # empty list for serach results
  # create an intermediate test list for atoms to be tested
  itl = dnts['close'] + dnts['dist']
  itl = my_clean_list(itl)
  if VerboseFlag==2:
    printf("Rings containing neighbours atoms\n")
    printf("    atom  smallest ring\n")
  for ne in itl:
    if VerboseFlag==2:
      printf("  %6s", MyLabel(ne))
      if len(ne.rings)>0:
        printf("  %i\n", len(ne.rings[0]))
      else:
        printf("  -\n")
    if len(ne.rings)>0 and len(ne.rings[0])<=10 : # Is the atom member of a ring?
      for at in ne.rings[0].atoms: # focus on the smallest ring
        if at not in qm:
          new_atoms.append(at)
  new_atoms = my_clean_list(new_atoms)
  qm += new_atoms
  qm  = my_clean_list(qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# distant neighbors                                                            #
# input   max          maximum distance in Angs                                #
#         fac          scaling factor for the sum of VdW radii                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def distant_neighbors(max, vdw_fac):
  # are the necessary global variables defined?
  gvar_list = ["qm", "dnts", "VerboseFlag"]
  for gvar in gvar_list:
    if gvar not in globals():
      printf("Variable %s NOT in globals.\n", gvar)
      exit()
    elif VerboseFlag==2: printf("Variable %s in globals.\n", gvar)
  global qm, dnts
  # the actual search
  new_atoms = []  # empty list for serach results
  bnd_cnt = 0      # counter for new bonds
  vdw_sum = 0.0    # sum of the VdW radii of a long distance contact

  if VerboseFlag == 2:
    printf("Distant neighbour search\n")
    printf("  <= %3.1f A and <= %.0f%% of the VdW radius\n", max, vdw_fac*100)
  for a1 in dnts['cent']: # loop over all center atoms
    for a2 in mol.atoms: # loop over all atoms in the molecule
      if a1 == a2 : continue
      dist = atom_dist(a1, a2)
      vdw_sum =  (1.0 * a1.vdw_radius)  # The type of the vdw_property appears to be
      vdw_sum += (1.0 * a2.vdw_radius)  # context dependent. I force it to a float.
      vdw_sum *= vdw_fac                # Apply cut-off factor
      if dist < max and a2 not in qm:
        if VerboseFlag == 2:
          printf("  %s%i",    a1.atomic_symbol, a1.index)
          printf( "-%s%i",    a2.atomic_symbol, a2.index)
          printf("  %6.4f",   dist)
          printf("  %6.4f",   vdw_sum)
          if dist<=vdw_sum: printf("  Passed\n")
          else: printf("  too long\n")
        if dist<=vdw_sum:
          new_atoms.append(a2)
          # Adding new bonds messes with the CDS assignment of rings. The new 
          # bonds can create small 3-member rings looping back to the metal 
          # atoms. To reactivate the option for new bonds change False to True
          if InMolBonds(mol, a1, a2)==-1 and bool(True):
            mol.add_bond(1, a1, a2)
            bnd_cnt += 1
  new_atoms = my_clean_list(new_atoms)
  qm += new_atoms
  qm = my_clean_list(qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# find close neighbors                                                         #
# input   none                                                                 #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
#         dnts         dictionary holding intermediate results                 #
################################################################################

def close_neighbors():
  # are the necessary global variables defined?
  gvar_list = ["qm", "dnts", "VerboseFlag"]
  for gvar in gvar_list:
    if gvar not in globals():
      printf("Variable %s NOT in globals.\n", gvar)
      exit()
    elif VerboseFlag==2: printf("Variable %s in globals.\n", gvar)
  global qm, dnts
  # the actual search
  new_atoms = []  # empty list for serach results
  for mc in dnts['cent']:
    for at in mc.neighbours:
      new_atoms.append(at)
      qm.append(at)
  new_atoms = my_clean_list(new_atoms)
  qm = my_clean_list(qm)
  return(new_atoms)

################################################################################
# functions to to help with the QM/MM separation                               #
# find metal center                                                            #
# input   at_num       atomic number of the lightest metal atom                #
# output  new_atoms    new atoms found and added to the QM core                #
# globals VerboseFlag  int to controll the verbosity of the output             #
#         qm           list with all atoms in the quantum core                 #
################################################################################

def metal_center(at_num):
  # are the necessary global variables defined?
  gvar_list = ["qm", "VerboseFlag"]
  for gvar in gvar_list:
    if gvar not in globals():
      printf("Variable %s NOT in globals.\n", gvar)
      exit()
    elif VerboseFlag==2: printf("Variable %s in globals.\n", gvar)
  global qm
  # the actual search
  new_atoms = []  # empty list for serach results
  for at in mol.atoms:
    if at.is_metal and at.atomic_number >= at_num:
      new_atoms.append(at)
  new_atoms = my_clean_list(new_atoms)
  qm += new_atoms
  qm = my_clean_list(qm)
  return(new_atoms)

################################################################################
# functions for the compartimensatiom of the  QM/MM separation                 #
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
  if len(flist)>0:
    printf("  Write %s\n", fname)
    list_xyz(flist, fname, mytxt)
  else:
    printf("  Skip  %s\n", fname)
  return bool(True)

################################################################################
# functions for the compartimensation of the  QM/MM separation                 #
# read input file                                                              #
# global VerboseFlag  controll the output to stdout                            #
# input  inp_name     name of the input file with the path to it               #
# output mol          CSD molecule object                                      #
#        obmol        OpenBabel molecule object                                #
################################################################################

def read_input_file(inp_name):

  # is VerboseFlag globally defined
  if "VerboseFlag" not in globals():
    print("Variable 'VerboseFlag' not globally defined")
    exit()

  # split input names into useful token
  root, ext = os.path.splitext(inp_name)
  path, file = os.path.split(inp_name)
  path += '/'
  base = root.replace(path, '')
  if VerboseFlag == 2:
    print("ext  ", ext)
    print("path ", path)
    print("base ", base)

  # check if the file exists
  if not os.path.isfile(inp_name):
    printf("file %s not found\n", inp_name)
    exit()

  # handle the input file by extention
  if VerboseFlag > 0: printf("Process input file %s\n", inp_name)
  input_type = ext.lower()
  if input_type not in ['.xyz', '.mol2']:
    printf("Wrong input file type\n")
    exit()

  # handel xyz file by converting them to mol2
  # If a mol2 file with same name is available, I will read the available 
  # file. If no substitute is available, the xyz file will be converted
  # into a mol2 file. Finally, the name of the input file will be updated
  # so that the input can be processed in the standard way.
  if input_type==".xyz":
    if VerboseFlag > 0: printf("  process xyz input file (depreciated)\n")
    mol2_inp_name = inp_name.replace(ext, ".mol2")
    if os.path.isfile(mol2_inp_name):
      if VerboseFlag > 0:
        printf("  %s (easier to use file format) found\n", mol2_inp_name)
      # update input information
      inp_name = mol2_inp_name
      ext = ".mol2"
      input_type = ".mol2"
    else:
      if VerboseFlag > 0:
        printf("  use Python module OpenBabel to convert %s%s\n", base, ext)
      obConversion = ob.OBConversion()
      obConversion.SetInAndOutFormats("xyz", "mol2")
      obmol = ob.OBMol()
      obConversion.ReadFile(obmol, inp_name)
      if VerboseFlag > 0:
        printf("  write converted file %s\n", mol2_inp_name)
      obConversion.WriteFile(obmol, mol2_inp_name)
      if VerboseFlag > 0:
        printf("  switch to %s\n", mol2_inp_name)
      if VerboseFlag == 2:
        printf("  delete ObenBabel object (rebuild later)\n")
      del obConversion, obmol
      # update input information
      inp_name = mol2_inp_name
      ext = ".mol2"
      input_type = ".mol2"

  # process the preferred mol2 input file
  # This part is in a conditional block so that the number of input
  # file types can be extended later.
  if input_type==".mol2":
    if VerboseFlag > 0:
      printf("  process mol2 input file (preferred)\n")
      printf("  read input file directly into a CSD object\n")
    mol_reader = io.MoleculeReader(inp_name)
    mol = mol_reader[0]
    if VerboseFlag > 0:
      printf("  standarize bonds\n")
    mol.assign_bond_types()
    mol.standardise_aromatic_bonds()
    mol.standardise_delocalised_bonds()

  # At this point I should have a CSD molecule object.
  # build auxiliary OpenBabel molecule object by reading the mol2 file
  if VerboseFlag > 0:
    printf("  build supplementary OpenBabel object from mol2 file\n")
  obConversion = ob.OBConversion()
  obConversion.SetInAndOutFormats("mol2", "mol2")
  obmol = ob.OBMol()
  obConversion.ReadFile(obmol, inp_name)
  if VerboseFlag > 0:
    printf("  compare the molecular geometries of both objects\n")
  csdanz = len(mol.atoms)
  obanz  = obmol.NumAtoms()
  if VerboseFlag==2:
    printf("  Number of atoms: CSD %i   OB %i\n", csdanz, obanz)
  if csdanz != obanz:
    print("The number of atoms in both molecule objects (CSD:  %i, OB: %i) don't match",
          csdanz, obanz)
    exit()
  else:
    if VerboseFlag>0:
      printf("    the number of atoms %i in both objects match\n", csdanz)
  if VerboseFlag>0:
    printf("    compare geometries line by line ...\n")
  for n in range(csdanz):
    obatom = obmol.GetAtom(n+1)
    csdnum = mol.atoms[n].atomic_number
    obnum  = obatom.GetAtomicNum()
    csdX = mol.atoms[n].coordinates.x
    csdY = mol.atoms[n].coordinates.y
    csdZ = mol.atoms[n].coordinates.z
    obX  = obatom.GetX()
    obY  = obatom.GetY()
    obZ  = obatom.GetZ()
    deltaX = abs(csdX-obX)
    deltaY = abs(csdY-obY)
    deltaZ = abs(csdZ-obZ)
    if VerboseFlag==2:
      printf("    %3i", n)
      printf(  "  %3i", csdnum)
      printf(  "  %3i", obnum)
      printf(  "  %10.6f", deltaX)
      printf(  "  %10.6f", deltaY)
      printf(  "  %10.6f", deltaZ)
      printf("\n")
    if csdnum != obnum:
      printf("The atomic numbers (CSD %i, OB %i) for entry #%i don't match\n",
            csdnum, obnum, n)
      exit()
    if deltaX >= 0.0001:
      printf("The X coordinates (CSD %.6f, OB %.6f) for entry #%i don't match\n",
            csdX, obX, n)
      exit()
    if deltaY >= 0.0001:
      printf("The Y coordinates (CSD %.6f, OB %.6f) for entry #%i don't match\n",
            csdY, obY, n)
      exit()
    if deltaZ >= 0.0001:
      printf("The Z coordinates (CSD %.6f, OB %.6f) for entry #%i don't match\n",
            csdZ, obZ, n)
      exit()
  if VerboseFlag>0: printf("    the geometries of both objects match\n")
  # return both molecule objects
  return(mol, obmol)

################################################################################
# functions for the compartimensation of the  QM/MM separation                 #
# read input file                                                              #
# global VerboseFlag  controll the output to stdout                            #
# input  mol          CSD mol                                                  #
#        smarts_list  list with SMARTS strings for the functionl groups        #
# output func_groups  list of lists containing the atoms of detected groups    #
################################################################################

def find_functional_groups(mol, smarts_list):
  func_groups = []
  # is VerboseFlag globally defined
  if "VerboseFlag" not in globals():
    print("Variable 'VerboseFlag' not globally defined")
    exit()
  # loop over list of smarts for functional groups
  if VerboseFlag>0:
    printf("Search for functional groups\n")
  for smarts in smarts_list:
    # create a substructure form the smarts string
    substrcuct = ccdc.search.SMARTSSubstructure(smarts)
    # search for the substructure
    substructure_search = SubstructureSearch()
    substructure_search.add_substructure(substrcuct)
    matches = substructure_search.search(mol)
    # check and report matches
    if matches:
      if VerboseFlag>0:
        printf("  Found %i %s group(s) in the molecule.\n", len(matches), smarts)
      cnt=0
      for match in matches:
        if VerboseFlag>0: printf("    %2i", cnt)
        fgroup=[]
        for at in match.match_atoms():
          str = "%s-%i" % (at.label, at.index)
          if VerboseFlag>0: printf("  %-6s", str)
          fgroup.append(at)
        if VerboseFlag>0: printf("\n")
        cnt+=1
        func_groups.append(fgroup)
    else:
      if VerboseFlag>0: printf("  No %s group(s) found.\n", smarts)
  if VerboseFlag>0:
    printf("  %i functional group(s) found.\n", len(func_groups))
  return(func_groups)

################################################################################
################################################################################
##                                                                            ##
##       MM   MM   AAAAA   II   NN   N   CCCCC   OOOOO   DDD     EEEE         ##
##       M M M M   A   A   II   N N  N   C       O   O   D  D    E            ##
##       M  M  M   AAAAA   II   N  N N   C       O   O   D   D   EEE          ##
##       M     M   A   A   II   N    N   C       O   O   D  D    E            ##
##       M     M   A   A   II   N   NN   CCCCC   OOOOO   DDD     EEEE         ##
##                                                                            ##
################################################################################
################################################################################


################################################################################
# parse comand line parameter (this should do for now)                         #
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

clean_glob_var(['my_arg', 'my_arglist', 'argpos', 'argcnt'])

################################################################################
# clean up old files                                                           #
################################################################################

if VerboseFlag > 0: printf("Clean up old files\n")
for file in glob.glob("stp*xyz"):
  os.remove(file)
if os.path.exists("./combi.xyz"): os.remove("./combi.xyz")

################################################################################
# split the command line argument for the input file                           #
################################################################################

mol, obmol = read_input_file(inp_name)

################################################################################
# identify functional groups                                                   #
################################################################################


# list of smart strings for the search of functional groups
# later rewrite this code to read the smarts form file
smarts_list=[]
smarts_list.append("[CX2]#N")                        # nitrile
smarts_list.append("[CX3](=O)[O]")                   # carboxylate
smarts_list.append("[CX3](=O)[H]")                   # aldehyde
smarts_list.append("[cX3][F,Cl,Br,I]")               # aryl halide
smarts_list.append("[#6X3](~O)~[#6X3]~[#6X3](~O)")   # enolate
smarts_list.append("[#6X3](=O)[#7X3](~H)")           # peptide
smarts_list.append("[#7X2]~[#7X2]~[#7X1]")           # azide

# create lookup list for functional goups in the molecule
func_groups = find_functional_groups(mol, smarts_list)

################################################################################
# the actual separation code                                                   #
################################################################################

if VerboseFlag > 0:
  printf("Start the actual QM/MM partioning process\n")

# find metal centers
stp_cnt = 1
dnts['cent'] = metal_center(19) # start with potassium
dnts.update()
summarize_step(stp_cnt, dnts['cent'], qm, "metal center")

# find the direct neighbours to the metal center
stp_cnt += 1
dnts['close'] = close_neighbors()
dnts.update()
summarize_step(stp_cnt, dnts['close'], qm, "direct neighbours")

# additional search based on distance, use the VdW radii for validity
stp_cnt += 1
dnts['dist'] = distant_neighbors(1.0, 0.8)
dnts.update()
summarize_step(stp_cnt, dnts['dist'], qm, "distance based neighbours")

# add the atoms of rings containing direct neighbours
stp_cnt += 1
dnts['nrings'] = neighbor_rings()
dnts.update()
summarize_step(stp_cnt, dnts['nrings'], qm, "rings containing neighbours atoms")

# first round of aromatic ring joined to inner rings
stp_cnt   += 1    # increase QM file counter
dnts['cArings'] = inner_aromatic_rings()
dnts.update()
summarize_step(stp_cnt, dnts['cArings'], qm, "inner aromatic rings")

printf("\nEntering search loop\n")
dnts['conju']   = []
dnts['gArings'] = []
dnts.update()
hc = []  # help/intermediate variable for conjugated chains
ha = []  # help/intermediate variable for aromatic rings
cnt = 0
while bool(True):
  cnt+=1
  num_new_at = 0
  # Test step to check for conjugated chanins
  stp_cnt += 1
  hc = conjugated_chains()
  num_new_at += len(hc)
  dnts['conju'] = my_clean_list(dnts['conju'] + hc)
  dnts.update()
  summarize_step(stp_cnt, hc, qm, "conjugated chains")
  # grow the number of joined aromatic ring
  stp_cnt   += 1
  ha = grow_aromatic_rings()
  num_new_at += len(ha)
  dnts['gArings'] += ha
  dnts.update()
  summarize_step(stp_cnt, ha, qm, "iterative search for aromatic rings")
  if num_new_at==0:
    printf("Round %i:  no new atoms\n", cnt)
    break
  else:
    printf("Round %i:  %2i new atoms\n\n", cnt, num_new_at)
printf("Leaving search loop\n\n")
del hc, ha, cnt

# looking for 1 atom links between QM atoms
stp_cnt   += 1
dnts['aLink'] = one_atom_links()
dnts.update()
summarize_step(stp_cnt, dnts['aLink'], qm, "1 atom links")

# looking for 2 atom links between the rings
stp_cnt   += 1
dnts['aaLink'] = two_atom_links()
dnts.update()
summarize_step(stp_cnt, dnts['aaLink'], qm, "2 atom links")

# adding single hydrogen atoms
stp_cnt   += 1
dnts['hatom'] = lone_H_atoms()
dnts.update()
summarize_step(stp_cnt, dnts['hatom'], qm, "H atoms")

# create mm layer
stp_cnt   += 1
mm = []
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
