"""
Author:	Nikko Stewart <nicho4@umbc.edu>
Python implementation of BEND by Dickerson & Goodsell
		Bending and curvature calculations in B-DNA.
		Goodsell DS; Dickerson RE; Nucleic Acids Res 22; 5497-503 (1994)
"""
import math
from numpy import *
sequencefile = open('sequence.dat', 'r')
profilelog = open('profile.log', 'w')
anglesfile = open('angles.dat', 'r')

sequence = []
isequence = []
"""Angles"""
twist = array([ (float(0), float(0), float(0), float(0)), 
						(float(0), float(0), float(0), float(0)), 
						(float(0), float(0), float(0), float(0)), 
						(float(0), float(0), float(0), float(0)) 
					])
roll = array([ 	(float(0), float(0), float(0), float(0)), 
						(float(0), float(0), float(0), float(0)), 
						(float(0), float(0), float(0), float(0)), 
						(float(0), float(0), float(0), float(0)) 
				])
tilt = array([ 	(float(0), float(0), float(0), float(0)), 
						(float(0), float(0), float(0), float(0)), 
						(float(0), float(0), float(0), float(0)), 
						(float(0), float(0), float(0), float(0)) 
				])

twistsum = float(0)
dx = 0
dy = 0
x = []
y = []
xave = []
yave = []
#bend = []
bend = zeros(500)
bendscale = 0
curve = []
curvescale = 0
pi = float(3.14159)
sc = 3.0

"""Read in the sequence from the file 'sequence.dat' """
def readNucleotides(sequence):
	n = 0
	"""No header lines, so read in everything"""
	for line in sequencefile:
		"""Read in entire line except '\n'"""
		for c in line:
			if (c != '\n'):
				sequence.append(c)
				n += 1
	print 'Number of bases: %d' % (n)
	return n	#Return number of sequence read in

"""Convert the sequence from the sequence to integers"""
def nucleotideToInt(isequence):
	i = 0
	
	"""BEND used A=1, B=2, .... but only because Fortran begins indexing at 1"""
	for nuc in sequence:
		if (nuc == 'A'):
			isequence[i] = 0
		elif (nuc == 'T'):
			isequence[i] = 1
		elif (nuc == 'G'):
			isequence[i] = 2
		elif(nuc == 'C'):
			isequence[i] = 3
		i += 1
	
	return
	
def readTwist(twist):
	"""4 lines of angles, with 4 angles per line"""
	pi = float(3.14159)
	for i in range(0, 4):
		line = anglesfile.readline()
		inputs = [float(z) for z in line.split(',')]
		for j in range(0, 4):
			twist[i][j] = float(inputs[j])
			#print "twist[%d][%d] = %.1f" % (i, j, twist[i][j])
			
	"""Convert each twist value to radians"""
	for i in range(0, 4):
		for j in range(0, 4):
			temp = float(twist[i][j])
			twist[i][j] = math.radians(temp)
			#print "twist[%d][%d] = %f" % (i, j, twist[i][j])

	return
"""Read in the roll angles"""
def readRoll(roll):
	"""4 lines of angles, with 4 angles per line"""
	for i in range(0, 4):
		"""Comma is delimiter"""
		temp = anglesfile.readline().split(',')
		"""Get rid of useless decimals and the newline"""	
		temp = [t.replace('\n', '') for t in temp]
		
		"""Print for debugging"""
		#print temp[0], temp[1] , temp[2], temp[3]
		"""Put each column of the line in angles"""
		for j in range(0, 4):
			roll[i][j] = float(temp[j])
	return

def readTilt(tilt):
	"""4 lines of angles, with 4 angles per line"""
	for i in range(0, 4):
		"""Comma is delimiter"""
		temp = anglesfile.readline().split(',')
		"""Get rid of useless decimals and the newline"""	
		temp = [t.replace('\n', '') for t in temp]

		"""Put each column of the line in angles"""
		for j in range(0, 4):
			tilt[i][j] = temp[j]
	return

for i in range(0, 500):
	#twist.insert(i, 0)
	#bend.append(float(0))
	x.insert(i, 0)
	y.insert(i, 0)
	xave.insert(i, 0)
	yave.insert(i, 0)
	curve.insert(i, 0)

"""Read in the sequence and convert the sequence to integers"""
n = readNucleotides(sequence)
"""Print sequence as formatted in file"""
s = ''
i = 0
while (i < n):
	for j in range(0, 10):
		s += sequence[i + j]
	print s
	s = ''
	i += 10
for i in range(0, n):
	isequence.insert(i, 0)
nucleotideToInt(isequence)
"""Ignore first 3 lines of the file"""
for i in range(0, 3):
	anglesfile.readline()

"""Read the input files and store the angles in matrices"""
readTwist(twist)
readRoll(roll)
readTilt(tilt)

temp = anglesfile.readline().split(',')
rbend, rcurve, bendscale, curvescale = (float(s) for s in temp)

if (rbend > 5):
	rbend = 5

twistsum = float(0)
"""Calculate trajectory of helix axis"""
for i in range(0, n-1):
	twistsum += twist[isequence[i]][isequence[i+1]]

	#print 'twistsum: %5f' % (twistsum)
	dx = (float(roll[isequence[i]][isequence[i+1]]) * math.sin(twistsum)) +  \
			(float(tilt[isequence[i]][isequence[i+1]]) * math.sin(twistsum - (pi/2)))
	dy = (float(roll[isequence[i]][isequence[i+1]]) * math.cos(twistsum)) +  \
			(float(tilt[isequence[i]][isequence[i+1]]) * math.cos(twistsum - (pi/2)))
	x.insert(i+1, x[i] + dx)
	y.insert(i+1, y[i] + dy)
	
"""Calculate average of trajectory over 10 base pairs"""
for i in range(5, n-5):
	rxsum = float(0);
	rysum = float(0);
	
	for j in range(-4, 5):
		rxsum += x[i+j]
		rysum += y[i+j]
	
	rxsum += (x[i+5]/float(2))
	rysum += (y[i+5]/float(2))
	rxsum += (x[i-5]/float(2))
	rysum += (y[i-5]/float(2))
	
	xave.insert(i, rxsum/float(10))
	yave.insert(i, rysum/float(10))
	

for i in range(int(rbend), int(n - rbend)):
	"""bend.insert(i, float(math.sqrt(	\
						(x[int(i+rbend)] - x[int(i-rbend)]) * (x[int(i+rbend)] - x[int(i-rbend)]) + \
						(y[int(i+rbend)] - y[int(i-rbend)]) * (y[int(i+rbend)] - y[int(i-rbend)])	\
						) * bendscale))"""
	bend[i] = (float(math.sqrt(	\
						(x[int(i+rbend)] - x[int(i-rbend)]) * (x[int(i+rbend)] - x[int(i-rbend)]) + \
						(y[int(i+rbend)] - y[int(i-rbend)]) * (y[int(i+rbend)] - y[int(i-rbend)])	\
						) * bendscale))

for i in range(int(rcurve+5), int(n-rcurve-6)):
	curve.insert(i, math.sqrt(	\
						(xave[i + int(rcurve)] - xave[i - int(rcurve)]) *	\
						(xave[i + int(rcurve)] - xave[i - int(rcurve)]) +	\
						(yave[i + int(rcurve)] - yave[i - int(rcurve)]) * 	\
						(yave[i + int(rcurve)] - yave[i - int(rcurve)])		\
						));
i = 0;
s =''
"""Print out the elements in each list, 50/line"""
while ((i + 50) < n):
	for j in range(i, i + 50):
		s += str(int(curve[j]/10))
	print s
	s = ''
	for j in range(i, i + 50):
		s += str(((int(curve[j] - (10 * int(curve[j] / 10))))))
	print s
	s = ''
	for j in range(i, i+50):
		s += sequence[j]
	print s, '\n'
	s = ''
	i += 50

"""Print out remaining nucleotides"""
remain = n - i;
for j in range(i, i + 50):
	s += str((int(curve[j])/10))
print s
s = ''
for j in range(i, i + 50):
	s += str((int(curve[j] - (10 * int(curve[j] / 10)))))
print s
s = ''
for j in range(i, i+remain):
	s += sequence[j]
print s, '\n'