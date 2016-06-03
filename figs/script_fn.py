import glob
import re
import os

files_to_process = glob.glob("conf*.png")

re_pattern = r"conf(\d+)(r?)\.png"

for filename in files_to_process:
    res = re.search( re_pattern, filename )
    if res.group(2):
        new_filename = "conf%03d.png" % (int(res.group(1))*2+1,)
    else:
        new_filename = "conf%03d.png" % (int(res.group(1))*2,)
    print "%s --> %s" % (filename, new_filename)
    os.rename(filename, new_filename)


