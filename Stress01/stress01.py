# 19.09.2024 - Building a data base for the stress tests
# - Starting with ccdc_t01, but use my experience from ccdc_q13.py

# Load modules for Python coding
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import math                     # mathematical  functions

# Load modules and prepare to read the CCS data base
import ccdc.search

from ccdc          import io                    # read crystal structures
from ccdc.search   import Search
from ccdc.search   import TextNumericSearch     # search by text or numbers
from ccdc.search   import SubstructureSearch    # search for substructures

from ccdc.molecule import Molecule              # Build a molecule
from ccdc.molecule import Atom                  # Atomic data and properties
from ccdc.molecule import Bond                  # Bond properties

# Fake printf and fprintf as used in C
def printf(format, *args):
  sys.stdout.write(format % args)
def fprintf(stream, format, *args):
  stream.write(format % args)

################################################################################
# Set up environment for the data base search                                  #
################################################################################

basedir = os.path.dirname(os.path.realpath(__file__))
resdir = basedir + "/Results"
printf("%-16s: ", "Base directory")
printf("%s\n",  basedir)
printf("%-16s: ", "Res directory")
printf("%s\n",  resdir)

if os.path.isdir(resdir):
  print("Directory for results exists")
else:
  print("Directory for results does NOT exist -> Create")
  os.makedirs(resdir)

################################################################################
# Search the CCDC database                                                     #
################################################################################

# Define target transition metal
target="Fe"
# list of unwanted transtion metals
tm=["Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y" , "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg"]
# remove Fe from the list of transition metals
tm.remove(target)
csd_reader = io.EntryReader('CSD')
cnt=0
for entry in csd_reader:
  if not entry.is_organometallic: continue
  mol = entry.molecule
  # Skipping all unwanted entries
  if len(mol.atoms) < 50: continue
  if (target+"1") not in mol.formula: continue
  if any(x in mol.formula for x in tm): continue
  # Analyze survivors
  cnt+=1
  printf("%7i  ", cnt)
  printf("%-10s  ", entry.identifier)
  printf("%4i  ", len(mol.atoms))
  printf("%s  ", mol.formula)
  printf("\n")
  if cnt==10: break
exit()


# Initialize CSD database connection
csd_reader = io.EntryReader('CSD')

# Search for entries with Fe atom and more than 100 atoms
search = Search.Search()
search.add_element('Fe')  # Iron atom

# Set filter for sixfold coordination
six_fold_coordination = Search.Parameters.CoordinationNumber(6)
search.add_coordinated_by('Fe', six_fold_coordination)

# Filter by number of atoms (greater than 100 atoms)
entries_with_iron_sixfold = []
for entry in csd_reader:
    if entry.molecule.atom_count > 100:  # More than 100 atoms
        # Find iron atoms and check if sixfold coordinated
        for atom in entry.molecule.atoms:
            if atom.atomic_symbol == 'Fe' and len(atom.neighbours) == 6:
                entries_with_iron_sixfold.append(entry.identifier)
                break  # No need to check further once condition is met

# Output the filtered entries
print(f'Found {len(entries_with_iron_sixfold)} entries with more than 100 atoms and sixfold coordinated iron:')
for entry_id in entries_with_iron_sixfold:
    print(entry_id)


exit()


# build a structure search by name
text_numeric_search = TextNumericSearch()
text_numeric_search.add_compound_name('aspirin')
identifiers = [h.identifier for h in text_numeric_search.search()]
identifiers = sorted(set(identifiers))
# print(identifiers)

# not all entries have a melting point
# printf("%i Enties found\n", len(identifiers))
# for i in range(len(identifiers)):
#     e = csd_reader.entry(identifiers[i])
#     printf("%3i  %-10s  %s\n", i, identifiers[i], e.melting_point)

cnt = 0
for ident in identifiers:
  e = csd_reader.entry(ident)
  mol = csd_reader.molecule(ident)
  if e.melting_point and hasattr(mol.atoms[0].coordinates, "x"):
    cnt += 1
    fname = "entry%02i.mol2" % cnt
    printf("%12s  %-8s  %s\n", fname, e.identifier, e.chemical_name)
    with io.MoleculeWriter(fname) as mol_writer: mol_writer.write(mol)
    fname = "entry%02i.xyz" % cnt
    out = open(fname, "w")
    fprintf(out, "%s\n", len(mol.atoms))
    fprintf(out, "%s\n", e.chemical_name)
    for n in range(len(mol.atoms)):
      if hasattr(mol.atoms[n].coordinates, "x"):
        fprintf(out, "%-2s  ", mol.atoms[n].atomic_symbol)
        fprintf(out, "%10.6f  ", mol.atoms[n].coordinates.x)
        fprintf(out, "%10.6f  ", mol.atoms[n].coordinates.y)
        fprintf(out, "%10.6f\n", mol.atoms[n].coordinates.z)
    out.close()

    
    
exit()