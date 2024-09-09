# 08.08.2024 - Getting started
# Read include funcionality to handle disorder
# The identifier is CIMHIT
# https://www.ccdc.cam.ac.uk/structures/Search?Doi=10.1038%2Fs41586-023-05821-2&DatabaseToSearch=Published

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
else:
  printf("  No disorder\n")

exit()
