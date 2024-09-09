# 07.08.2024 - Getting started
# Read one file (ferrocene) and analyze it

# Load modules for Python coding
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system

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

# load tools to handle CCSD database
from ccdc import io                       # read crystal structures
csd_reader = io.EntryReader('CSD')
from ccdc.search import TextNumericSearch # search by text or numbers

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