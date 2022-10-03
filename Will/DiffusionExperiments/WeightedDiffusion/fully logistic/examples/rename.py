import os

for root, dirs, files in os.walk('.'):
    for name in files:
        if '.png' in name:
            if len(name) == 22:
                newname = name[:15] + '0' + name[15:]
                print(name, newname)
                os.rename(name, newname)
            if len(name) == 24:
                newname = name[:15] + name[16:]
                print(name, newname)
                os.rename(name, newname)