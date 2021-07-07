import serial
import numpy as np
import keyboard
import sys

np.set_printoptions(threshold=sys.maxsize)
arduino = serial.Serial(port = 'COM5', baudrate = 38422, timeout = .1)

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
    filename = f'{value}.csv'
    file1 = open(filename, "w")

    value = input("Press Enter to start recording\n")
    print("Press Space to stop recording")
    
    while True:
        data = read_serial().strip()
        if data:
            strdata = str(data)
            file1.writelines(strdata[2:len(strdata) - 1])
            file1.writelines("\n")
        if keyboard.is_pressed('space'):
            print(f'Stopped recording and closed file {filename}')
            file1.close()
            break




