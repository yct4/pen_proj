import serial
import numpy as np
import keyboard
import sys
import matlab.engine as matlab
import re

np.set_printoptions(threshold=sys.maxsize)
arduino = serial.Serial(port = 'COM5', baudrate = 38422, timeout = .1)
eng = matlab.start_matlab()

def read_serial():
    data = arduino.readline()
    return data

# init loop
print("Start initialization")
for i in range(10):
    data = read_serial().strip()
    if data:
        print(data)
print("End initialization")

# main loop
while True:
    value = input("Enter name of output file:\n")
    print(f'You entered {value}')
    filename = f'data\\{value}.csv'
    file1 = open(filename, "w")

    temp = input("Press Enter to start recording\n")
    print("Press Space to stop recording")
    
    while True:
        data = read_serial().strip()
        if data:
            strdata = str(data)
            file1.writelines(strdata[2:len(strdata) - 1])
            file1.writelines("\n")
        if keyboard.is_pressed('space'):
            print(f'Stopped recording and closed file {filename}. File is stored in ./data')
            file1.close()
            break

    # py2matlab
    print("running matlab to get chaincode")
    chaincodes = np.array(eng.quat9_0705(value,nargout=1)).astype(int)
    chaincodes_str = re.sub('[\[\]]', '', np.array2string(chaincodes, separator=','))
    print("printing chaincodes")
    print(chaincodes_str)
    chaincodes_len = len(chaincodes)
    print(f'\nchaincodes # rows = {chaincodes_len}')
    
    # manually classify chaincodes
    actletters = ""
    j = chaincodes_len
    for i in range(chaincodes_len):
        actletters += input(f'Enter actual classification of letter {j}:\n').strip() + '\n'
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



