import numpy as np
import keyboard
import sys
import matlab.engine as matlab

np.set_printoptions(threshold=sys.maxsize)

# init loop
print("Start initialization")
value = input("Enter name of output file:\n")
print(f'You entered {value}')
eng = matlab.start_matlab()
eng.quat9_0705(value,nargout=0)
print('iteration 2')
eng.quat9_0705(value,nargout=0)
print("End initialization")

while True:
    i = 0

# main loop
#while True:
#    value = input("Enter name of output file:\n")
#    print(f'You entered {value}')
#    filename = f'{value}.csv'
#    file1 = open(filename, "w")
#
#    #value = input("Press Enter to start recording\n")
#    #print("Press Space to stop recording")
#    file1.writelines("this is a test string\n")
#    file1.close()





