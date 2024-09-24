################################################################################
# 24.09.2024
# - This Puthon code holds all the functions to be imported
################################################################################

# Import basic Python modules
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import math                     # mathematical  functions

# Import CCDC modules (CSD Python API)
from ccdc          import io                    # read crystal structures
csd_reader = io.EntryReader('CSD')

import ccdc.search
from ccdc.search   import TextNumericSearch     # search by text or numbers
from ccdc.search   import SubstructureSearch    # search for substructures

from ccdc.molecule import Molecule              # Build a molecule
from ccdc.molecule import Atom                  # Atomic data and properties
from ccdc.molecule import Bond                  # Bond properties

# Import OpenBabel to assist with atom classification
from openbabel import openbabel as ob           # basic OpenBabel

# Import my modules
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q15_CCDC')
import globals as gb  # short for globals

################################################################################
# general purpose function                                                     #
# fake printf and fprintf as used in C                                         #
################################################################################

def printf(format, *args):
  sys.stdout.write(format % args)
def fprintf(stream, format, *args):
  stream.write(format % args)

################################################################################
# functions for the compartimensation of the  QM/MM separation                 #
# read input file                                                              #
# global VerboseFlag  controll the output to stdout                            #
# input  inp_name     name of the input file with the path to it               #
# output mol          CSD molecule object                                      #
#        obmol        OpenBabel molecule object                                #
################################################################################

def read_input_file(inp_name):

  # are the global variable initialized
  if "gb" not in globals():
    print("Global variables not initialized")
    exit()

  # split input names into useful token
  root, ext = os.path.splitext(inp_name)
  path, file = os.path.split(inp_name)
  path += '/'
  base = root.replace(path, '')
  if gb.VerboseFlag == 2:
    print("ext  ", ext)
    print("path ", path)
    print("base ", base)

  # check if the file exists
  if not os.path.isfile(inp_name):
    printf("file %s not found\n", inp_name)
    exit()

  # handle the input file by extention
  if gb.VerboseFlag > 0: printf("Process input file %s\n", inp_name)
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
    if gb.VerboseFlag > 0: printf("  process xyz input file (depreciated)\n")
    mol2_inp_name = inp_name.replace(ext, ".mol2")
    if os.path.isfile(mol2_inp_name):
      if gb.VerboseFlag > 0:
        printf("  %s (easier to use file format) found\n", mol2_inp_name)
      # update input information
      inp_name = mol2_inp_name
      ext = ".mol2"
      input_type = ".mol2"
    else:
      if gb.VerboseFlag > 0:
        printf("  use Python module OpenBabel to convert %s%s\n", base, ext)
      obConversion = ob.OBConversion()
      obConversion.SetInAndOutFormats("xyz", "mol2")
      obmol = ob.OBMol()
      obConversion.ReadFile(obmol, inp_name)
      if gb.VerboseFlag > 0:
        printf("  write converted file %s\n", mol2_inp_name)
      obConversion.WriteFile(obmol, mol2_inp_name)
      if gb.VerboseFlag > 0:
        printf("  switch to %s\n", mol2_inp_name)
      if gb.VerboseFlag == 2:
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
    if gb.VerboseFlag > 0:
      printf("  process mol2 input file (preferred)\n")
      printf("  read input file directly into a CSD object\n")
    mol_reader = io.MoleculeReader(inp_name)
    mol = mol_reader[0]
    if gb.VerboseFlag > 0:
      printf("  standarize bonds\n")
    mol.assign_bond_types()
    mol.standardise_aromatic_bonds()
    mol.standardise_delocalised_bonds()

  # At this point I should have a CSD molecule object.
  # build auxiliary OpenBabel molecule object by reading the mol2 file
  if gb.VerboseFlag > 0:
    printf("  build supplementary OpenBabel object from mol2 file\n")
  obConversion = ob.OBConversion()
  obConversion.SetInAndOutFormats("mol2", "mol2")
  obmol = ob.OBMol()
  obConversion.ReadFile(obmol, inp_name)
  if gb.VerboseFlag > 0:
    printf("  compare the molecular geometries of both objects\n")
  csdanz = len(mol.atoms)
  obanz  = obmol.NumAtoms()
  if gb.VerboseFlag==2:
    printf("  Number of atoms: CSD %i   OB %i\n", csdanz, obanz)
  if csdanz != obanz:
    print("The number of atoms in both molecule objects (CSD:  %i, OB: %i) don't match",
          csdanz, obanz)
    exit()
  else:
    if gb.VerboseFlag>0:
      printf("    the number of atoms %i in both objects match\n", csdanz)
  if gb.VerboseFlag>0:
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
    if gb.VerboseFlag==2:
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
  if gb.VerboseFlag>0: printf("    the geometries of both objects match\n")
  # return both molecule objects
  return(mol, obmol)

