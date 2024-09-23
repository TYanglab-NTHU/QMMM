import sys
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q14_CCDC/WebExmaple')
import settings

def stuff():
    settings.myList.append('hey')

def more_stuff():
    print(settings.myList)

def extra_stuff():
    settings.myList.append("Murkel")