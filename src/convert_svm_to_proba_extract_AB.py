#!/usr/bin/env python3

import sys
import math

#def sigmoid(x):
#  return 1 / (1 + math.exp(-x))


### Input parameters

# out = array of SVM outputs
# target = array of booleans: is ith example a positive example
# prior1 = number of positive examples
# prior2 = number of negative examples

### Output parameteres
# A, B = parameters of sigmoid


# Define variable


prior0 = 0.0
prior1 = 0.0


# Load SVM output
svm_output = []
svm_output_file = open(sys.argv[1], "r")
while 1:
    line = svm_output_file.readline()
    if line == "":
        break
    line = line.rstrip()
    svm_score = float(line)
    svm_output.append(svm_score)
svm_output_file.close()


# Load SVM target
target = []
gene_site_list = []
svm_target_file = open(sys.argv[2], "r")
while 1:
    line = svm_target_file.readline()
    if line == "":
        break
    tab = line.split()
    if tab[0] != "#":
        svm_target = tab[0]
        if svm_target[0] == "+":
            prior0 = prior0+1
            target.append(1)
        elif svm_target[0] == "-":
            prior1 = prior1+1
            target.append(0)
        #print tab
        gene_site = tab[-1][1:]
        #print gene_site
        gene_site_list.append(gene_site)
svm_target_file.close()

#sys.exit()

# Main
A = 0.0
B = math.log((prior0+1)/(prior1+1))
hiTarget = (prior1+1)/(prior1+2)
loTarget = 1/(prior0+2)
lambda_value = 1e-3
olderr = 1e300


##print hiTarget
##print loTarget
##print lambda_value
##print olderr


# pp = temp array to store current estimate of probability of
# examples all pp array elements to (prior1+1)/(prior0+prior1+2)
pp = []
for i in range(len(target)):
    pp.append((prior1+1)/(prior0+prior1+2))
#print pp


count = 0
for it in range(100):
    # Initiliase values
    a, b, c, d, e = 0, 0, 0, 0, 0

    # First, compute Hessian and gradient of error function
    # with respect to A and B.
    for i in range(len(target)):
        if target[i]:
            t = hiTarget
        else:
            t = loTarget

        d1 = pp[i]-t
        d2 = pp[i]*(1-pp[i])

        a = a+(svm_output[i]*svm_output[i]*d2)
        b = b+d2
        c = c+svm_output[i]*d2
        d = d+svm_output[i]*d1
        e = e+d1
    # If gradient is really tiny, then stop
    if (abs(d) < 1e-9) and (abs(e) < 1e-9):
        break
    oldA = A
    oldB = B
    err = 0.0

    # Loop until goodness of fit increase
    while 1:
        det = (a+lambda_value)*(b+lambda_value)-c*c
        if det == 0:  # if determinant of Hessian is zero
            lambda_value = lambda_value*10 # increase stabiliser
            continue
        A = oldA + ((b+lambda_value)*d-c*e)/det
        B = oldB + ((a+lambda_value)*e-c*d)/det

        # Now, compute the goodness of fit
        err = 0.0
        for i in range(len(target)):
            p = 1/(1+math.exp(svm_output[i]*A+B))
            pp[i] = p

            # At this step, make sur log(0) returns = -200
            err = err-t*math.log(p)+(1-t)*math.log(1-p)
        if err < olderr*(1+1e-7):
            lambda_value = lambda_value*0.1
            break
        # Error did not decrease, increase stabiliser by factor of 10
        # and try again
        lambda_value = lambda_value*10
        if lambda_value > 1e6: # Something is broken. Give up
            break

    diff = float(err)-float(olderr)
    scale = 0.5*(err+olderr+1)
    if (diff > -1e-3*scale) and (diff < 1e7*scale):
        count = count+1
    else:
        count = 0
    olderr = err
    if count == 3:
        break

print(A, B)

# Problem with Value domain math error
# => Use Sigmoid instead

# sys.exit()
output = open(sys.argv[3], "w")
output.write(str(A)+"\t"+str(B)+"\n")
output.close()
