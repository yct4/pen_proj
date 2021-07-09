import numpy as np
import keyboard
import sys
import re
from sklearn import svm
import os

np.set_printoptions(threshold=sys.maxsize)

directory = r'C:\Users\nico\pen_proj\chaincodes'

# read from chaincodes files
chaincodes_train = []
classes_train = []
count = 0
directory_train = directory + r'\train'
for filename in os.listdir(directory_train):
    if filename.endswith(".csv"):
        print(os.path.join(directory, filename))
    else:
        continue

#print(directory_train)
for filename in os.listdir(directory_train):
    if filename.endswith("chaincodes.csv"):
        chaincodes = np.genfromtxt(directory_train + '\\' + filename, delimiter=',', dtype='int', usecols=range(29))
        if (count < 1):
            chaincodes_train = chaincodes
        else:
            chaincodes_train = np.concatenate([chaincodes_train, chaincodes])
    elif filename.endswith("classify.csv"):
        classes = np.genfromtxt(directory_train + '\\' + filename, delimiter=',', dtype='str')
        if (count < 1):
            classes_train = classes
        else:
            classes_train = np.concatenate([classes_train, classes])
    count = count + 1

classes_a = classes_train != '!'
chaincodes_train = chaincodes_train[classes_a]
classes_train = classes_train[classes_a]

print("chaincodes_train: ")
print(chaincodes_train)
print("classes_train: ")
print(classes_train)


# read from chaincodes files
chaincodes_test = []
classes_test = []
count = 0
directory_test = directory + r'\test'
for filename in os.listdir(directory_test):
    if filename.endswith(".csv"):
        print(os.path.join(directory, filename))
    else:
        continue

#print(directory_test)
for filename in os.listdir(directory_test):
    if filename.endswith("chaincodes.csv"):
        chaincodes = np.genfromtxt(directory_test + '\\' + filename, delimiter=',', dtype='int', usecols=range(29))
        if (count < 1):
            chaincodes_test = chaincodes
        else:
            chaincodes_test = np.concatenate([chaincodes_test, chaincodes])
    elif filename.endswith("classify.csv"):
        classes = np.genfromtxt(directory_test + '\\' + filename, delimiter=',', dtype='str')
        if (count < 1):
            classes_test = classes
        else:
            classes_test = np.concatenate([classes_test, classes])
    count = count + 1

classes_a = classes_test != '!'
chaincodes_test = chaincodes_test[classes_a]
classes_test = classes_test[classes_a]

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



