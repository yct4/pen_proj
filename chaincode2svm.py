import numpy as np
import keyboard
import sys
import re
from sklearn import svm

np.set_printoptions(threshold=sys.maxsize)

# init loop
#print("Start initialization")
#value = input("Enter name of input file (must be in directory ./chaincodes):\n")
#print(f'You entered {value}')

train_data = ['test1', 'letter_test_0705']
test_data = ['letter_test_0705_2']

# read from chaincodes files
chaincodes_train = []
classes_train = []
for value in train_data:
    filename = f'chaincodes\\train\{value}_chaincodes.csv'
    filename_classes = f'chaincodes\\train\{value}_chaincodes_classify.csv'
    chaincodes = np.genfromtxt(filename, delimiter=',', dtype='int', usecols=range(29))
    classes = np.genfromtxt(filename_classes, delimiter=',', dtype='str')

    classes_a = classes != '!'
    chaincodes = chaincodes[classes_a]
    classes = classes[classes_a]

chaincodes_train = chaincodes
classes_train = classes

chaincodes_test = chaincodes_train

# SVM rbf kernel
#rbf_svc = svm.SVC(kernel='rbf')
clf = svm.SVC(kernel='rbf', decision_function_shape='ovr')
fit_svm = clf.fit(chaincodes_train, classes_train)
dec = clf.decision_function(chaincodes_train)
dec.shape[1] # 4 classes: 4*3/2 = 6

pred = clf.predict(chaincodes_test)

#print(clf)
#print(fit_svm)
#print(dec)
#print(dec.shape[1])
#print(pred)

#clf.decision_function_shape = "ovr"
#dec = clf.decision_function([[1]])
#dec.shape[1] # 4 classes

#chaincodes = np.array(eng.quat9_0705(value,nargout=1)).astype(int)
##chaincodes_str = np.array2string(chaincodes_resub, separator=',')
#chaincodes_str = re.sub('[\[\]]', '', np.array2string(chaincodes, separator=','))
#print("printing chaincodes")
#print(chaincodes_str)
#
#chaincodes_len = len(chaincodes)
#print(f'\nchaincodes # rows = {chaincodes_len}')
#
#actletters = ""
#for i in range(chaincodes_len):
#    actletters += '\n' + input(f'Enter actual classification of letter {i+1}:\n')
#actletters_len = len(actletters) 
#print(f'\nActual Letter Classifications:\n{actletters[0:actletters_len]}\n')
#
## WRITE TO CHAINCODES FILE
#filename = f'chaincodes\{value}_chaincodes.csv'
#file1 = open(filename, "w")
#file1.writelines(chaincodes_str)
#file1.writelines(actletters[0:actletters_len])
#file1.close()
#
#print("Chaincodes file: {value}_chaincodes.csv is saved in ./chaincodes")
#print("End initialization")

#while True:
#    i = 0

