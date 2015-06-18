import time
import sys

for i in range(11):
    time.sleep(1)
    sys.stdout.write("\r|"+"-"*i+"."*(10-i)+"|%d%%" % (i*10))
    sys.stdout.flush()