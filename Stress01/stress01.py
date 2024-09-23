# 23.09.2024
# - picking random files without repetions

import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import random                   # create random number
import time                     # gather timing info

sys.path.insert(0, '/dicos_ui_home/tlankau/QMMM/Q14_CCDC')
from qmmm_func import printf
from qmmm_func import read_input_file
from qmmm_func import metal_center

# define vital variables
file_dir = "/dicos_ui_home/tlankau/QMMM/Query02/QRes-134592/Results"
stress_anz = 1

def mol2_num_atoms(file_name):
  name=""
  anz=0
  f = open(file_name,"r")
  linelist = f.readlines()
  cnt=0
  for line in linelist:
    if "@<TRIPOS>MOLECULE" in line : break
    cnt += 1
  name = linelist[cnt+1].strip()
  anz = linelist[cnt+2].strip()
  anz = anz.split()
  anz = int(anz[0])
  return(name, anz)  

done = []
total_time = 0
for cnt in range(1, stress_anz+1, 1):
  printf("%3i  ", cnt)
  while bool(True):
    file = random.choice(glob.glob(file_dir+'/*.mol2'))
    if file not in done: break
  done.append(file)
  name, numat = mol2_num_atoms(file)
  printf("%-8s  %3i  ", name, numat)
  start_time = time.time()

  ##############################################################################
  # Do the QMM/MM separation                                                   #
  ##############################################################################
  

  
  
  
  ##############################################################################
  # Leave the QMM/MM separation                                                #
  ##############################################################################
   
  end_time   = time.time()
  time_diff = end_time - start_time
  total_time += time_diff
  printf("%.1f", time_diff)
  printf("\n")
printf("Avergae time per QM/MM separation %.1f\n", total_time/stress_anz)  