# SIMPLE PYTHON SCRIPT FOR LINE PLOTS

import matplotlib.pyplot as plt
import sys
import numpy as np

f = open("{}".format(sys.argv[1]),"r")
ft = f.readlines()
f.close()

if len(sys.argv) == 5:
	# YXX PLOT
	col1 = int(sys.argv[2]) - 1
	col2 = int(sys.argv[3]) - 1
	col3 = int(sys.argv[4]) - 1

	k = 0
	ft1 = ft[k].split()
	xl = ft1[col1].strip("#")
	yl1 = ft1[col2]
	yl2 = ft1[col2].strip("\n")
	k = k + 1

	list1 = []
	list2 = []
	list3 = []
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = float(ft1[col1].strip("\n"))
		t1 = round(t1,2)
		t2 = float(ft1[col2].strip("\n"))
		t2 = round(t2,2)
		t3 = float(ft1[col3].strip("\n"))
		t3 = round(t3,2)
		list1.append(t1)
		list2.append(t2)
		list3.append(t3)
		k = k + 1

	x = np.array(list1)
	y1 = np.array(list2)
	y1max = np.amax(y1)
	y1min = np.amin(y1)
	y2 = np.array(list3)
	y2max = np.amax(y2)
	y2min = np.amin(y2)

	if y1min < y2min:
		min1 = y1min - 2.0
	else:
		min1 = y2min - 2.0

	if y1max > y2max:
		max1 = y1max + 2.0
	else:
		max1 = y2max + 2.0

	fig, ax1 = plt.subplots()

	color = 'tab:red'
	ax1.set_ylabel('{}'.format(xl), fontsize=10, fontweight='bold')
	ax1.set_xlabel('{}'.format(yl1), color=color, fontsize=10, fontweight='bold')
	ax1.plot(y1, x, color=color)
	ax1.set_xlim([min1,max1])
	ax1.tick_params(axis='x', labelcolor=color)

	ax2 = ax1.twiny()  

	color = 'tab:blue'
	ax2.set_xlabel('{}'.format(yl2), color=color, fontsize=10, fontweight='bold')  
	ax2.plot(y2, x, color=color)
	ax2.tick_params(axis='x', labelcolor=color)
	ax2.set_xlim([min1,max1])
	fig.tight_layout() 
	plt.show()

elif len(sys.argv) == 4:

	col1 = int(sys.argv[2]) - 1
	col2 = int(sys.argv[3]) - 1

	k = 0
	ft1 = ft[k].split()
	xl = ft1[col1].strip("#")
	yl1 = ft1[col2]
	k = k + 1

	list1 = []
	list2 = []
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = float(ft1[col1].strip("\n"))
		t1 = round(t1,2)
		t2 = float(ft1[col2].strip("\n"))
		t2 = round(t2,2)
		list1.append(t1)
		list2.append(t2)
		k = k + 1

	x = np.array(list1)
	y1 = np.array(list2)

	fig, ax1 = plt.subplots()

	ax1.set_ylabel('{}'.format(yl1), fontsize=10, fontweight='bold')
	ax1.set_xlabel('{}'.format(xl), fontsize=10, fontweight='bold')
	ax1.scatter(x, y1, color='black')
	ax1.plot(x, y1, color='black')
	
	fig.tight_layout() 
	plt.show()

	
		
