# THIS CODE HELPS IN PLOTING THE OUTPUT CONTRAST MAPS FROM AMBAT_AT

#••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# THIS CODE REQUIRES PYTHON 3 WITH SCIPY AND NUMPY LIBRARIES
# PLOTING IS DONE BY:
		# CONVERTING POLAR OUTPUT FROM AMBAT TO CARTESIAN COORDINATES
		# Delaunay Triangulation
		# LinearNDInterpolation :: TO GET THE CONTINOUS HEAT MAP FROM DISCRETE AMBAT INPUT
		# COUNTOUR PLOTING

# NOTE :: VARIABLE NAME npoints CONTROL THE SMOOTHNESS OF THE GRAPH, DEFAULT IS 50, WHICH RESULTS IN 2500 GRID POINTS

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

# CODE SYNTAX

# python PLOTS_AMBAT_AT.py $filename $column(of property) 3D(optional for 3D protein profiles)

# •••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••

from pylab import *
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator as IN
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
from scipy.interpolate import griddata

def frange(start,stop,step):
	i = start
	while i < stop:
		yield(i)
		i += step

def heat_map():

	npoints = 50

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	col = int(sys.argv[2])

	# CHANGING THE COORDINATES TO CARTESIAN COORDINATES

	k = 1
	list1 = []
	list2 = []
	while k < len(ft):
		ft1 = ft[k].split()
		radius = float(ft1[0])
		angle = float(ft1[1])
		z = float(ft1[(col-1)].strip("\n"))

		x = radius * (math.cos((angle*pi/180.0)))
		y = radius * (math.sin((angle*pi/180.0)))

		list3 = [x,y]
		list1.append(list3)

		list2.append(z)
		
		k = k + 1

	points = np.array(list1)
	values = np.array(list2)
	tri = Delaunay(points)
	inter = IN(tri,values,fill_value=-1.00)

	# GETTING THE RANGE OF R AND THETA

	ft1 = ft[1].split()
	rin = float(ft1[0].strip("\n")) 
	tin = float(ft1[1].strip("\n"))
	tin = 2*np.pi*tin/180.0
	ft1 = ft[(len(ft)-1)].split()
	rfi = float(ft1[0].strip("\n"))
	tfi = float(ft1[1].strip("\n"))
	tfi = 2*np.pi*tfi/180.0

	rstep = (rfi - rin) / npoints
	rstep = round(rstep,2)

	tstep = (tfi - tin) / npoints
	tstep = round(tstep,2)

	listx = []
	listy = []
	list2 = []
	list3 = []
	list4 = []
	for i in frange(rin,rfi,rstep):
		listx.append(i)
	for j in frange(tin,tfi,tstep):
		listy.append(j)

	for i in frange(rin,rfi,rstep):
		for j in frange(tin,tfi,tstep):
			X = i * (math.cos(j))
			Y = i * (math.sin(j))
			Z = float(inter(X,Y))
			list3.append(Z)

	values = np.array(list3)
	values.shape = (len(listx),len(listy))
	xvalue = np.array(listx)
	yvalue = np.array(listy)

	fig, axs = plt.subplots(1, subplot_kw=dict(projection='polar'))
	p2 = axs.contourf(yvalue,xvalue, values, cmap='RdBu')
	cbar = plt.colorbar(p2, ax=axs)
	label_cb = sys.argv[1].split(".")
	cbar.ax.set_ylabel('{}'.format(label_cb[0]), fontsize=10, fontweight='bold')
	#plt.ylim(10,30)
	plt.show()

def heat_map_3D():

	npoints = 50

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	col = int(sys.argv[2])

	# CHANGING THE COORDINATES TO CARTESIAN COORDINATES

	ft1 = ft[1].split()
	initial_z = ft1[2]
	final_z = initial_z

	# GETTING THE NUMBER OF SLICES
	k = 1
	count = 1
	ztemp =initial_z
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[2]
		if t1 != ztemp:
			ztemp = t1
			count = count + 1
		k = k + 1

	print("THERE ARE {} SLICES IN THE FILE. WHICH ONE YOU WANT TO PLOT?".format(count))
	slices = input()

	#rows = count - 2
	#fig, axs = plt.subplots(rows,2, subplot_kw=dict(projection='polar'))
	fig, axs = plt.subplots(1, subplot_kw=dict(projection='polar'))
	#nrows = rows
	#ncols = 2
	#x1 = 0
	#y1 = 0

	# GETTING THE RANGE OF R AND THETA

	ft1 = ft[1].split()
	rin = float(ft1[0].strip("\n")) 
	tin = float(ft1[1].strip("\n"))
	tin = 2*np.pi*tin/180.0
	ft1 = ft[(len(ft)-1)].split()
	rfi = float(ft1[0].strip("\n"))
	tfi = float(ft1[1].strip("\n"))
	tfi = 2*np.pi*tfi/180.0

	rstep = (rfi - rin) / npoints
	rstep = round(rstep,2)

	tstep = (tfi - tin) / npoints
	tstep = round(tstep,2)

	k = 1
	list1 = []
	list2 = []
	listz = []
	s = 1
	while k < len(ft):
		ft1 = ft[k].split()
		radius = float(ft1[0])
		angle = float(ft1[1])
		z = float(ft1[(col-1)].strip("\n"))
		final_z = ft1[2]
		if final_z == initial_z:
			x = radius * (math.cos((angle*pi/180.0)))
			y = radius * (math.sin((angle*pi/180.0)))

			list3 = [x,y]
			list1.append(list3)

			list2.append(z)
			intital_z = final_z
		else:
			points = np.array(list1)
			values = np.array(list2)
			tri = Delaunay(points)
			inter = IN(tri,values,fill_value=-1.00)

			listx = []
			listy = []
			list2 = []
			list3 = []
			list4 = []
			for i in frange(rin,rfi,rstep):
				listx.append(i)
			for j in frange(tin,tfi,tstep):
				listy.append(j)

			for i in frange(rin,rfi,rstep):
				for j in frange(tin,tfi,tstep):
					X = i * (math.cos(j))
					Y = i * (math.sin(j))
					Z = float(inter(X,Y))
					list3.append(Z)

			if s == int(slices) - 1:
				from_z = initial_z

			if int(slices) == 1:
				from_z = "lower_leaflet"

			if s == int(slices):
				values = np.array(list3)
				values.shape = (len(listx),len(listy))
				xvalue = np.array(listx)
				yvalue = np.array(listy)
				#p2 = axs[x1, y1].contourf(yvalue,xvalue, values, cmap='RdBu')
				#axs[x1, y1].set_title('Z = {}'.format(initial_z))
				p2 = axs.contourf(yvalue,xvalue, values, cmap='RdBu')
				axs.set_title('Z = {} to {}'.format(from_z,initial_z),loc="left",fontsize=10, fontweight='bold',pad=20.0)

				#cbar = plt.colorbar(p2, ax=axs[x1, y1])

			#if x1 < (nrows-1):
				#x1 = x1 + 1
			#else:
				#y1 = y1 + 1
	
			initial_z = final_z
			s = s + 1
			k = k - 1
			list1 = []
			list2 = []
			listz = []

		k = k + 1

	# FOR THE FINAL LEAFLET

	if s == int(slices):
		points = np.array(list1)
		values = np.array(list2)
		tri = Delaunay(points)
		inter = IN(tri,values,fill_value=-1.00)

		listx = []
		listy = []
		list2 = []
		list3 = []
		list4 = []
		for i in frange(rin,rfi,rstep):
			listx.append(i)
		for j in frange(tin,tfi,tstep):
			listy.append(j)

		for i in frange(rin,rfi,rstep):
			for j in frange(tin,tfi,tstep):
				X = i * (math.cos(j))
				Y = i * (math.sin(j))
				Z = float(inter(X,Y))
				list3.append(Z)

		if s == int(slices) - 1:
			from_z = initial_z

		if int(slices) == 1:
			from_z = "lower_leaflet"


		values = np.array(list3)
		values.shape = (len(listx),len(listy))
		xvalue = np.array(listx)
		yvalue = np.array(listy)
		#p2 = axs[x1, y1].contourf(yvalue,xvalue, values, cmap='RdBu')
		#axs[x1, y1].set_title('Z = {}'.format(initial_z))
		p2 = axs.contourf(yvalue,xvalue, values, cmap='RdBu')
		axs.set_title('Z = {} to {}'.format(from_z,initial_z),loc="left",fontsize=10, fontweight='bold',pad=20.0)

	
	cbar = plt.colorbar(p2, ax=axs)
	label_cb = sys.argv[1].split(".")
	cbar.ax.set_ylabel('{}'.format(label_cb[0]), fontsize=10, fontweight='bold')
	#plt.ylim(10,30)
	plt.show()


# CHECK IF INPUT IS COORECT

if len(sys.argv) < 3:
	print("ERROR IN INPUT")
	print("SYNTAX :: python PLOTS_AMBAT_AT.py $filename $column(of property) 3D(optional for 3D protein profiles)")
	quit()

try:
	arg1 = sys.argv[3]
except:
	arg1 = "1D"

if arg1 == "3D" or arg1 == "3d":
	heat_map_3D()
else:
	heat_map()






