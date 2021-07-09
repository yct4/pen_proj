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
count = 0
for value in train_data:
    filename = f'chaincodes\\train\{value}_chaincodes.csv'
    filename_classes = f'chaincodes\\train\{value}_chaincodes_classify.csv'
    chaincodes = np.genfromtxt(filename, delimiter=',', dtype='int', usecols=range(29))
    classes = np.genfromtxt(filename_classes, delimiter=',', dtype='str')

    classes_a = classes != '!'
    chaincodes = chaincodes[classes_a]
    classes = classes[classes_a]

    if (count < 1):
        chaincodes_train = chaincodes
        classes_train = classes
    else:
        chaincodes_train = np.concatenate([chaincodes_train, chaincodes])
        classes_train = np.concatenate([classes_train, classes])
    count = count + 1

print("chaincodes_train: ")
print(chaincodes_train)
print("classes_train: ")
print(classes_train)


chaincodes_test = []
classes_test = []
count = 0
for value in test_data:
    filename = f'chaincodes\\test\{value}_chaincodes.csv'
    filename_classes = f'chaincodes\\test\{value}_chaincodes_classify.csv'
    chaincodes = np.genfromtxt(filename, delimiter=',', dtype='int', usecols=range(29))
    classes = np.genfromtxt(filename_classes, delimiter=',', dtype='str')

    classes_a = classes != '!'
    chaincodes = chaincodes[classes_a]
    classes = classes[classes_a]

    if (count < 1):
        chaincodes_test = chaincodes
        classes_test = classes
    else:
        chaincodes_test = np.concatenate([chaincodes_test, chaincodes])
        classes_test = np.concatenate([classes_test, classes])
    count = count + 1

print("chaincodes_test: ")
print(chaincodes_test)
print("classes_test: ")
print(classes_test)

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
print(f'prediction: {pred}')
print(f'actual    : {classes_test}')



