import numpy as np
import keyboard
import sys
import matlab.engine as matlab
import re

np.set_printoptions(threshold=sys.maxsize)

# init loop
print("Start initialization")
value = input("Enter name of input file (must be in directory ./data):\n")
print(f'You entered {value}')
eng = matlab.start_matlab()
print("started matlab")
chaincodes = np.array(eng.quat9_0705(value,nargout=1)).astype(int)
#chaincodes_str = np.array2string(chaincodes_resub, separator=',')
chaincodes_str = re.sub('[\[\]]', '', np.array2string(chaincodes, separator=','))
print("printing chaincodes")
print(chaincodes_str)

chaincodes_len = len(chaincodes)
print(f'\nchaincodes # rows = {chaincodes_len}')

actletters = ""
j = chaincodes_len
for i in range(chaincodes_len):
    actletters = input(f'Enter actual classification of letter {j}:\n') + '\n' + actletters
    j -= 1
actletters_len = len(actletters) 
print(f'\nActual Letter Classifications:\n{actletters[0:actletters_len]}\n')

# WRITE TO CHAINCODES FILE
filename = f'chaincodes\{value}_chaincodes.csv'
file1 = open(filename, "w")
file1.writelines(chaincodes_str)
file1.close()

filename = f'chaincodes\{value}_chaincodes_classify.csv'
file1 = open(filename, "w")
file1.writelines(actletters[0:actletters_len])
file1.close()


print("Chaincodes file: {value}_chaincodes.csv is saved in ./chaincodes")
print("End initialization")

while True:
    i = 0

