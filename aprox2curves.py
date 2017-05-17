#! /usr/bin/python3
import os
import math
from scipy.interpolate import interp1d

N = 100					# number of points of approximation
f1 = 'data/exp/M186/p_p0.txt' 			# data for curve №1
f2 = 'data/theory/M186/upper.xy'					# data for curve №2
# x must be first column
n1 = 2						          # column for curve №1
n2 = 2						          # column for curve №2
scale1x = 1					        # scale for curve №1 along the x axis
scale2x = 1					        # scale for curve №2 along the x axis
scale1y = 638331.089701564	# scale for curve №1 along the y axis 638331.089701564
scale2y = 1					        # scale for curve №2 along the y axis

# function of extract of data
def extract_data(filename, column, scale):
    infile = open(filename, 'r')
    numbers = []
    for line in infile:
        words = line.split()
        if words[0] == '#':
            infile.readline()
        else:
            number = float(words[column])*scale
            numbers.append(number)
    infile.close()
    return numbers

# find X1 and X2
_f1 = open(f1, 'r')
_f2 = open(f2, 'r')
for line in _f1:
    words = line.split()
    if words[0] == '#':
        _f1.readline()
    else:
        x11 = [float(w) for w in _f1.readline().split()]
        x12 = [float(w) for w in _f1.readlines()[-1].split()]
for line in _f2:
    words = line.split()
    if words[0] == '#':
        _f2.readline()
    else:
        x21 = [float(w) for w in _f2.readline().split()]
        x22 = [float(w) for w in _f2.readlines()[-1].split()]
_f1.close()
_f2.close()
if x11[0]*scale1x > x21[0]*scale2x:
	x1 = x11[0]*scale1x
else:
	x1 = x21[0]*scale2x
if x12[0]*scale1x > x22[0]*scale2x:
	x2 = x22[0]*scale2x
else:
	x2 = x12[0]*scale1x

i = 0
lx = []
for _ in range (N):
	if i+1 < N:
		lx.append(x1+(x2-x1)/(N-1)*i)
	else:
		lx.append(x2)
	i+=1
_lx1 = extract_data(f1, 0, scale1x)
_lx2 = extract_data(f2, 0, scale2x)
_ly1 = extract_data(f1, n1-1, scale1y)
_ly2 = extract_data(f2, n2-1, scale2y)
ly1_f = interp1d(_lx1, _ly1)
#ly2_f = interp1d(_lx2, _ly2)
ly1 = ly1_f(lx)
#ly2 = ly2_f(lx)
z = 0
ly2 = []
for _ in range (N):
	if z < N/2:
		ly2.append(_ly2[0])
	else:
		ly2.append(_ly2[3])
	z+=1
curv = open('Curves.txt', 'w')
j = 0
for _ in range (N):
    curv.write(str(lx[j]) + '\t')
    curv.write(str(ly1[j]) + '\t')
    curv.write(str(ly2[j]) + '\n')
    j+=1
curv.close()

k = 0
dy = []
for _ in range (N):
    dy_ = (ly1[k]-ly2[k])
    dy.append(dy_)
    k+=1
dys = sum(dy)/N
l = 0
sig_i = []
otn_i = []
for _ in range (N):
    sig_i_=(dy[l]-dys)**2
    sig_i.append(sig_i_)
    otn_i_ = math.fabs(dy[l]/ly1[l])*100
    otn_i.append(otn_i_)
    l+=1
sigma = (sum(sig_i)/N)**0.5
otn = sum(otn_i)/N
print('Standard deviation: ', sigma)
print('Medium relative precision: ', otn, '%')
