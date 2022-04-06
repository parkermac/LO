"""
Code to test the indexing logic of extract_box_monthly.py.

"""
import numpy as np

## Start of month loop

#NFN = 8641
NFN = 8761

print('NFN = ' + str(NFN))

a = np.arange(NFN)
b = np.array_split(a, 12)
    
for c in b:
    print('%d:%d (%d)' % (c[0], c[-1],len(c)))