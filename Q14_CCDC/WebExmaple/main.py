import sys
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q14_CCDC/WebExmaple')
import settings
import sub_file

settings.init()           # Call only once
sub_file.stuff()          # Do stuff with global var
print(settings.myList[0]) # Check the result
settings.myList.append('Glubb')
sub_file.more_stuff()     # Do stuff with global var
sub_file.extra_stuff()    # Do stuff with global var
sub_file.more_stuff()     # Do stuff with global var
print(settings.myList)