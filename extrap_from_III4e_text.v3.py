
#from os import listdir
import sys
import math
import os
import glob
import getopt
import xml.dom.minidom
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LineString
from shapely import geometry


## Integrate the steps from findI4ePoints.py into this program
##extrapolate off of I4e points
##	-- mkdir 'deg_intersect_points', 'fid_point'
##	-- ls *extrapolated.txt | xargs -I {} bash -c "python3 /Users/nicoletatro/git/3dGVF/programs/findI4ePoints.py -i {}"
##		-- ls *.txt if no open isopters
##	-- python3 /Users/nicoletatro/git/3dGVF/programs/extrap_from_I4e_text.py 

def rotate_points(avgRotation,avgRRatio,mX,mY,xlist,ylist):
	xlist_rot = []
	ylist_rot = []

	for i in range(0,len(xlist)):
		cx = xlist[i] - mX
		cy = mY - ylist[i]
		ptheta = math.atan2((cy),cx)
		if (ptheta < 0):
			ptheta += 2*math.pi
		pr = math.hypot(cx,cy)
		gtheta = ptheta - avgRotation
		if (gtheta < 0):
			gtheta += 2*math.pi
		elif (gtheta > 2*math.pi):
			gtheta -= 2*math.pi
		gr = pr/avgRRatio
		xTest = avgRRatio*gr*math.cos(gtheta) + mX
		yTest = mY - avgRRatio*gr*math.sin(gtheta)
		xlist_rot.append(xTest)
		ylist_rot.append(yTest)
	return (xlist_rot, ylist_rot)

def read_in_points(annotation,avgRotation,avgRRatio,mX, mY):
	polyNode = annotation.getElementsByTagName('spline')[0]
	xNode = polyNode.getElementsByTagName('xpoints')[0]
	yNode = polyNode.getElementsByTagName('ypoints')[0]
	num_points = len(xNode.getElementsByTagName("int"))
	xValues = xNode.getElementsByTagName("int")
	yValues = yNode.getElementsByTagName("int")
	xlist = []
	ylist = []
	for i in range(num_points):
		x = int(xValues[i].firstChild.nodeValue)
		xlist.append(x)
		y = int(yValues[i].firstChild.nodeValue)
		ylist.append(y)
	xlist2, ylist2 = rotate_points(avgRotation,avgRRatio,mX, mY, xlist, ylist)
	for i in range(len(ylist2)):
		ylist2[i] = (ylist2[i] - mY)*(-1)
		ylist2[i] = ylist2[i] + mY
	return (xlist2, ylist2)

def read_in_points_ex(annotation,avgRotation,avgRRatio,mX, mY):
	polyNode = annotation.getElementsByTagName('spline')[0]
	xNode = polyNode.getElementsByTagName('xpoints')[0]
	yNode = polyNode.getElementsByTagName('ypoints')[0]
	num_points = len(xNode.getElementsByTagName("int"))
	xValues = xNode.getElementsByTagName("int")
	yValues = yNode.getElementsByTagName("int")
	xlist = []
	ylist = []
	for i in range(num_points):
		x = int(xValues[i].firstChild.nodeValue)
		xlist.append(x)
		y = int(yValues[i].firstChild.nodeValue)
		ylist.append(y)
	xlist2, ylist2 = rotate_points(avgRotation,avgRRatio,mX, mY, xlist, ylist)
	#for i in range(len(ylist2)):
	#       ylist2[i] = (ylist2[i] - mY)*(-1)
	#       ylist2[i] = ylist2[i] + mY
	return (xlist2, ylist2)

def collectIsopters(each,avgRotation,avgRRatio,mX,mY):
	annotations = each.getElementsByTagName("annotation")
	isopters = []
	isopters_open = []
	isopters_extrap = []

	numIII4e = 1
	numIII4e_open = 1
	numIII4e_extrap = 1

	#Isopters
	for annotation in annotations:
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'III4e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				isopters.append('III4e' + str(numIII4e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				numIII4e = numIII4e + 1
			if (name == 'III4e-extrap'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				isopters_extrap.append('III4e-extrap' + str(numIII4e_extrap))
				isopters_extrap.append(xlist2)
				isopters_extrap.append(ylist2)
				numIII4e_extrap = numIII4e_extrap + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'III4e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				isopters.append('III4e-open' + str(numIII4e_open))
				isopters.append(xlist2)
				isopters.append(ylist2)
				numIII4e_open = numIII4e_open + 1

	return(isopters, isopters_open, isopters_extrap)

def rotateFiducials(avgRotation,avgRRatio,mX,mY,tX,tY,rX,rY,bX,bY,lX,lY):
	fidX = [tX,rX,bX,lX]
	fidY = [tY,rY,bY,lY]
	fidX_pt = [90,0,270,180]
	fidY_pt = [70,90,70,90]
	rFidX = []
	rFidY = []
	rFidX_pt = []
	rFidY_pt = []
	fidComplete = []
	ptsComplete = []

	for i in range(0,len(fidX)):
		testX = fidX[i]
		testY = fidY[i]
		testX_pt = fidX_pt[i]
		testY_pt = fidY_pt[i]
		if (testX != 0 and testY != 0):
			rFidX.append(testX)
			rFidY.append(testY)
			rFidX_pt.append(testX_pt)
			rFidY_pt.append(testY_pt)
	for i in range(0,len(rFidX)):
		cx = rFidX[i] - mX
		cy = mY - rFidY[i]
		ptheta = math.atan2((cy),cx)
		if (ptheta < 0):
			ptheta += 2*math.pi
		pr = math.hypot(cx,cy)
		gtheta = ptheta - avgRotation
		if (gtheta < 0):
			gtheta += 2*math.pi
		elif (gtheta > 2*math.pi):
			gtheta -= 2*math.pi
		gr = pr/avgRRatio
		xTest = avgRRatio*gr*math.cos(gtheta) + mX
		yTest = mY - avgRRatio*gr*math.sin(gtheta)
	for i in range(len(rFidX)):
		fidComplete.append(rFidX[i])
		fidComplete.append(rFidY[i])
		ptsComplete.append(rFidX_pt[i])
		ptsComplete.append(rFidY_pt[i])
	return (fidComplete, ptsComplete)



def rotationCalculation(mX,mY,tX,tY,rX,rY,bX,bY,lX,lY):
	length = 0
	avgRotation = 0
	avgRRatio = 0
	dist_fid = [tX,tY,rX,rY,bX,bY,lX,lY]
	dist_loc = [90,70,0,90,270,70,180,90]
	fid = []
	loc = []

	for i in range(0,len(dist_fid)):
		test = dist_fid[i]
		point = dist_loc[i]
		if (test != 0):
			length = length + 1
			fid.append(test)
			loc.append(point)
	#Number of fiducial points annotated
	length = (length/2)
	#Set up conversion and rotation
	for i in range(0,len(fid),2):
		#Cartesian coordinates translation
		#print (fid[i], fid[i+1])
		cx = fid[i] - mX
		cy = mY - fid[i+1]
		#pixel polar coordinates
		ptheta = math.atan2((cy),cx)
		pr = math.hypot(cx,cy)
		#GVF coordinate conversion set up
		thetaDifference = (ptheta - (math.radians(loc[i])))
		while (abs(thetaDifference)>math.pi):
			if(thetaDifference<0):
				thetaDifference += 2*math.pi
			else:
				thetaDifference -= 2*math.pi
		## need to compare some epsilon less than 0 to account for variability.
		if(thetaDifference<-0.1):
				thetaDifference += 2*math.pi
		avgRotation += (thetaDifference/length)
		avgRRatio += (pr/loc[i+1])/length
	return (avgRotation, avgRRatio)




def fiducials_scale(each):
	bottomX = 0
	bottomY = 0
	topX = 0
	topY = 0
	leftX = 0
	leftY = 0
	rightX = 0
	rightY = 0
	middleX = 0
	middleY = 0
	num_fid = 0
	size = 0
	scale = 0
	distance = 0
	distance_alt = 0
	eye = ''

	annotations = each.getElementsByTagName("annotation")
	#Fiducial points
	for annotation in annotations:
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "90,70    "):
			topX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			topY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "270,70    "):
			bottomX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			bottomY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "180,90    "):
			leftX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			leftY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "0,90    "):
			rightX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			rightY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "0,0    "):
			middleX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			middleY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "Eye?"):
			eye = annotation.getElementsByTagName("value")[0].firstChild.nodeValue

	if (middleX == 0 and middleY == 0):
		if (num_fid == 4):
			middleX = (topX + bottomX)/2
			middleY = (leftY + rightY)/2
	if (middleX != 0 and middleY != 0):
		if (topX != 0 and topY != 0):
			distance = math.sqrt(((topX-middleX)*(topX-middleX))+((topY-middleY)*(topY-middleY)))
		elif (bottomX != 0 and bottomY != 0):
			distance = math.sqrt(((bottomX-middleX)*(bottomX-middleX))+((bottomY-middleY)*(bottomY-middleY)))
		elif (leftX != 0 and leftY != 0):
			distance = math.sqrt(((leftX-middleX)*(leftX-middleX))+((leftY-middleY)*(leftY-middleY)))
		elif (rightX != 0 and rightY != 0):
			distance = math.sqrt(((rightX-middleX)*(rightX-middleX))+((rightY-middleY)*(rightY-middleY)))
		else:
			distance = 0
	if (distance > 0):
		size = distance/300
		scale = 25*size

	return(scale,num_fid,middleX,middleY,topX,topY,rightX,rightY,bottomX,bottomY,leftX,leftY,eye)


def get_imagename(each):
	imagename = each.getElementsByTagName("imagename")
	oldFileName = imagename[0].firstChild.nodeValue
	#print("DEBUG: starting value of oldFileName = ", oldFileName)
	#newFileName = []
	oldFileName = oldFileName.split('/')
	if len(oldFileName) > 1:
		oldFileName = oldFileName[1]
	else:
		oldFileName = oldFileName[0]
	#newFileName = ''.join(oldFileName)
	#newFileName = newFileName.split('_')
	#if len(newFileName) > 1:
	#	newFileName = newFileName[1]
	#else:
	#	newFileName = newFileName[0]
	#newFileName = ''.join(newFileName)
	fileName = oldFileName
	#print("DEBUG: final value of fileName = ", fileName)
	return(fileName)


def parseAnnotation(inputfile):
	inputfile = inputfile
	dom = xml.dom.minidom.parse(inputfile)
	images = dom.getElementsByTagName("annotation-per-image")

	for each in images:
		#Name of the output PDF
		fileName = get_imagename(each)
		#print("DEBUG: get_imagename: %s --> %s" % (each, fileName))
		#sys.exit()
		#Check fiducials and find scaling factor
		scale, num_fid,mX,mY,tX,tY,rX,rY,bX,bY,lX,lY,eye = fiducials_scale(each)

		## remove this after testing!
		#return(mX, mY, scale)


		if (mX != 0 and mY != 0 and scale != 0):
			print (fileName)
			sufficientFid = "Yes"
			#Calculate average rotation and RRatio for
			avgRotation, avgRRatio = rotationCalculation(mX,mY,tX,tY,rX,rY,bX,bY,lX,lY)
			#Rotate fiducial points
			rotatedFids, pts = rotateFiducials(avgRotation,avgRRatio,mX,mY,tX,tY,rX,rY,bX,bY,lX,lY)
			#Will rotate points and flip y-axis
			isopters,isopters_open,isopters_extrap = collectIsopters(each,avgRotation,avgRRatio,mX,mY)
			#calc_points(fileName,scale,isopters,isopters_open,isopters_extrap,mX,mY,lX,lY,bX,bY,rX,rY,tX,tY,eye)
			## this is OK, but not great. OK because there should only ever be one image per XML file at this point in the process
			return(mX, mY, scale, isopters)


def readFIDFile(fid_file):
	#print("In readFIDFile: ", fid_file)
	f = open(fid_file, 'r')
	line1 = f.readline().rstrip("\n")
	line2 = f.readline().rstrip("\n")
	f.close()
	#print("LINE1: ", line1)
	#print("LINE2: ", line2)
	line1s1 = line1.split("\t")
	#print("LINE1S1[1]: ", line1s1[1])
	line1s2 = line1s1[1].split(",")
	mx = line1s2[0]
	my = line1s2[1]
	line2s1 = line2.split("\t")
	scale = line2s1[1]
	#print("MX = ", mx, "\nMY = ", my, "\nSCALE = ", scale)
	## re-set scale to be pixels per 1 radial degree
	scale = float(scale) / 5.833333
	#print("MX = ", mx, "\nMY = ", my, "\nSCALE = ", scale)
	return(mx, my, scale)

## process the intersect points to build the list of lines to add to the XML file
def readIntersectFile(mX, mY, ppd, isoIntercepts, nd):
	#print("In readIntersectFile: ", intersect_file)
	III4e_X = []
	III4e_Y = []
	V4e_X = []
	V4e_Y = []
	for i in range(0,360):
		theta = (270 - float(i) + 360) % 360
		#theta = i
		#print("II=", isoIntercepts[i][0])
		radial_dist = convertCartesianToPolar(isoIntercepts[i][0], isoIntercepts[i][1], float(mX), float(mY), ppd)
		#print("INPUT:", arr, "\tTHETA=", theta, "\tDIST=", radial_dist, "\n")
		#rel_I4e = radial_dist / nd['I4e'][theta]
		rel_III4e = radial_dist / nd['III4e'][theta]
		extrap_V4e = rel_III4e * nd['V4e'][theta]
		#print(rel_I4e, " --> ", extrap_V4e)
		#X1,Y1 = convertPolarToCartesian(theta, extrap_III4e, mX, mY, ppd, 0.0)
		X2,Y2 = convertPolarToCartesian(theta, extrap_V4e, mX, mY, ppd, 0.0)
		#III4e_X.append(X1)
		#III4e_Y.append(Y1)
		V4e_X.append(X2)
		V4e_Y.append(Y2)
	#print("III4e_X= ", III4e_X)
	#print("III4e_Y= ", III4e_Y)
	#print("V4e_X= ", V4e_X)
	#print("V4e_Y= ", V4e_Y)

	to_add = ""
	## Not needed for III4e-based extrapolation
	#to_add += "\t\t<annotation>\n\t\t\t<name>III4e-extrap-full                             </name>\n\t\t\t<value>III4e-extrap-full                             </value>\n"
	#to_add += "\t\t\t<controltype>SplineRoi</controltype>\n\t\t\t<color>\n\t\t\t\t<red>255</red>\n\t\t\t\t<green>0</green>\n\t\t\t\t<blue>0</blue>\n\t\t\t</color>\n"
	#to_add += "\t\t\t<probability>0</probability>\n\t\t\t<annotation-start>1495556692.619784</annotation-start>\n"
	#to_add += "\t\t\t<annotation-end>1495556895.923801</annotation-end>\n\t\t\t<roi>\n\t\t\t\t<spline>\n\t\t\t\t\t<npoints>360</npoints>\n"

	#to_add += "\t\t\t\t\t<xpoints>\n"
	#for i in III4e_X:
	#	to_add += "\t\t\t\t\t\t<int>" + str(round(i)) + "</int>\n"
	#to_add += "\t\t\t\t\t</xpoints>\n"
	#to_add += "\t\t\t\t\t<ypoints>\n"
	#for i in III4e_Y:
	#	to_add += "\t\t\t\t\t\t<int>" + str(round(i)) + "</int>\n"
	#to_add += "\t\t\t\t\t</ypoints>\n"
	#to_add += "\t\t\t\t</spline>\n\t\t\t</roi>\n\t\t</annotation>\n"

	#to_add = "\n"

	to_add += "\t\t<annotation>\n\t\t\t<name>V4e-extrap-full                             </name>\n\t\t\t<value>V4e-extrap-full                                </value>\n"
	to_add += "\t\t\t<controltype>SplineRoi</controltype>\n\t\t\t<color>\n\t\t\t\t<red>255</red>\n\t\t\t\t<green>0</green>\n\t\t\t\t<blue>0</blue>\n\t\t\t</color>\n"
	to_add += "\t\t\t<probability>0</probability>\n\t\t\t<annotation-start>1495556692.619784</annotation-start>\n"
	to_add += "\t\t\t<annotation-end>1495556895.923801</annotation-end>\n\t\t\t<roi>\n\t\t\t\t<spline>\n\t\t\t\t\t<npoints>360</npoints>\n"

	to_add += "\t\t\t\t\t<xpoints>\n"
	for i in V4e_X:
		to_add += "\t\t\t\t\t\t<int>" + str(round(i)) + "</int>\n"
	to_add += "\t\t\t\t\t</xpoints>\n"
	to_add += "\t\t\t\t\t<ypoints>\n"
	for i in V4e_Y:
		to_add += "\t\t\t\t\t\t<int>" + str(round(i)) + "</int>\n"
	to_add += "\t\t\t\t\t</ypoints>\n"
	to_add += "\t\t\t\t</spline>\n\t\t\t</roi>\n\t\t</annotation>\n"

	#print(to_add)

	return(to_add)


def convertPolarToCartesian(theta, radius, mX, mY, ppd, avgRotation):
	## convert the angle to radians
	radian = theta * 2.0 * math.pi / 360
	## add the rotation back into the set of points
	radian = radian + avgRotation
	X = math.cos(radian)
	Y = math.sin(radian)
	X = float(mX) + (X * radius * ppd);
	Y = float(mY) - (Y * radius * ppd);
	#print("THETA=", theta, "\tRADIAN=", radian, "\tX=", X, "\tY=", Y, "\n")
	return(X,Y)
	## These are going to be AFTER the rotation from the original image, so they will be rotated TWICE!

def convertCartesianToPolar(pX,pY, mX, mY, ppd):
	#print("In convertC2P")
	dist = math.sqrt( (pX-mX)**2 + (pY-mY)**2) / ppd
	return(dist)

def normativeData():
	nd_I4e = { 0: 84, 15: 82, 30: 78, 45: 69, 60: 61.5, 75: 57, 90: 55, 105: 56.5, 120: 57.5, 135: 59.5, 150: 60, 165: 60.5, 180: 60, 195: 60, 210: 60, 225: 60, 240: 60.5, 255: 62.5, 270: 68, 285: 72.5, 300: 74.5, 315: 77.5, 330: 80, 345: 83.5}
	nd = {'I4e': nd_I4e}
	nd_III4e = { 0: 87, 15: 84.5, 30: 81, 45: 72.5, 60: 64.5, 75: 60, 90: 58, 105: 59.5, 120: 60.5, 135: 62.5, 150: 63, 165: 63, 180: 63, 195: 62.5, 210: 62.5, 225: 62.5, 240: 63, 255: 65.5, 270: 71, 285: 74.5, 300: 77.5, 315: 80.5, 330: 83, 345: 86}
	nd['III4e'] = nd_III4e
	nd_V4e = { 0: 90, 15: 87, 30: 84, 45: 76, 60: 67.5, 75: 63, 90: 61.5, 105: 62.5, 120: 63.5, 135: 65, 150: 66, 165: 66, 180: 66, 195: 65, 210: 65, 225: 65, 240: 65.5, 255: 68, 270: 74, 285: 76.5, 300: 80.5, 315: 83.5, 330: 86, 345: 89}
	nd['V4e'] = nd_V4e
	return(nd)

def main(argv):
	## process arguments
	inputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:",["input="])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			usage()
			sys.exit()
		elif opt in ("-i", "--input"):
			inputfile = arg
			print ('Input file is ', inputfile)
	if(inputfile == ''):
		usage()
		sys.exit()

	print("DEBUG: After processing command-line arguments")

	## Load normative data
	normData = normativeData()
	#print(normData)

	print("DEBUG: After loading normative data")

	mX, mY, scale, isopters = parseAnnotation(inputfile)
	print("DEBUG: isopters =", isopters)
	ppd = float(scale) / 5.833333


### from calc_points in find...
	#isopter to polygon
	polygon = ''

	if (len(isopters) > 0):
		#point_list = []
		## if more than one isopter, we want to use the largest one
		max_range_x = 0
		max_range_y = 0
		max_isopter = 0
		if(len(isopters) > 3):
			num_isopters = int(len(isopters) / 3)
			print("NUM_ISO=", num_isopters)
			#for i in range(0,len(isopters),3):
			for ii in range(0,num_isopters):
				min_x = 10000
				max_x = -10000
				#min_y = 10000
				#max_y = -10000
				iso_idx = (ii * 3) + 1
				for j in range(0,len(isopters[iso_idx])):
					if isopters[iso_idx][j] > max_x:
						max_x = isopters[iso_idx][j]
					if isopters[iso_idx][j] < min_x:
						min_x = isopters[iso_idx][j]
					#if isopters[iso_idx+1][j] > max_y:
					#	max_y = isopters[iso_idx+1][j]
					#if isopters[iso_idx+1][j] < min_y:
					#	min_y = isopters[iso_idx+1][j]
				if (max_x - min_x) > max_range_x:
					max_range_x = max_x - min_x
					max_isopter = ii
		print("MAX_ISOPTER=", max_isopter)
		iso_idx = (max_isopter * 3) + 1
		point_list = []
		#print (isopters[i])
		for j in range(0,len(isopters[iso_idx])):
			point_list.append((isopters[iso_idx][j], isopters[iso_idx+1][j]))
			#print((isopters[i+1][j], isopters[i+2][j]))
		polygon = Polygon(point_list)






	## Building a "complete" normalized dataset.
	## this should really be done once, and hard-coded into the normData 
	new_normData = {}
	for each in normData:
		#print (each)
		new_normData[each]={}
		## creating a local copy of normData for the isopter specified in "each"
		degrees = []
		distances = []
		for degree in normData[each]:
			degrees.append(degree)
			distances.append(normData[each][degree])
		#print (degrees)
		#print (distances)
		## linear interpolation of normalized data between the points in normData
		### goes from every 15 degrees in normData to every degree in normal_360_*
		normal_360_deg = []
		normal_360_dist = []
		## this handles 0 - 344
		for i in range(0,len(degrees)-1):
			deg_diff = degrees[i+1] - degrees[i]
			dist_diff = distances[i+1] - distances[i]
			deg_increment = 1
			dist_increment = dist_diff/deg_diff
			for j in range(0,deg_diff):
				new_deg = degrees[i] + j* deg_increment
				new_dist = distances[i] + dist_increment*j
				#print (new_deg, new_dist)
				normal_360_deg.append(new_deg)
				normal_360_dist.append(new_dist)
		## this handles 345 - 359
		deg_increment = ((normData[each][0] - normData[each][345])/15)
		for i in range(1, 16):
			new_deg = 344 + i
			new_dist = normData[each][345] + deg_increment*i
			#print (new_deg, new_dist)
			normal_360_deg.append(new_deg)
			normal_360_dist.append(new_dist)
		for i in range(0,len(normal_360_deg)):
			new_normData[each][normal_360_deg[i]] = normal_360_dist[i]


	## Process the file, extrapolating the I4e isopter
	file = inputfile
	## determine output filename
	arr = file.split(".")
	#print (arr)
	base = arr[0]
	updated_file_name = base + ".updated.txt"
	## read through the inputfile, and determine the lines to add for V4e and III4e based upon the I4e
	## - this logic is from readIntersectFile

	## open output file
	print("Opening output file:", updated_file_name)
	f2 = open(updated_file_name, "w")
	## read-through and re-write all of the lines prior to the end of the annotations
	remainder = ""
	remainder_flag = 0
	with open(file) as file_in:
		for line in file_in:
			if remainder_flag == 0:
				if not ('</annotations>' in line):
					line.rstrip("\n")
					# TODD print (line)
					f2.write(line)
				else:
					## read in the remainder of the file into xml_remainder, to print at the end
					remainder_flag = 1
					remainder = line

			else:
				remainder += line



## from find_I4e... 
	## can I use ppd instead of chunk? That would save needing to pass back bX and bY
	## I THINK it's safe to not use chunk
	#d_x = bX - mX
	#d_y = bY - mY
	#deg_length = math.sqrt((d_x*d_x)+(d_y*d_y))
	#chunk = float(deg_length/7)
	chunk = ppd * 10
	### TODD - I think chunk needs to be negative
	if(chunk > 0):
		chunk = chunk * -1.0


	## this is for the f2 that was the intersection points
	#f2 = open(pts_dir+fileName+"_pts.txt","w")

	## Need to load the polygon object with the data for the I4e...





	## After we have made it through the for loop, isoIntercepts should have entries for every degree, with both X and Y components
	isoIntercepts = {}
	for i in range(0,360):
		isoIntercepts[i] = {}
		x1 = mX
		y1 = mY
		x2 = (12*chunk*(math.sin(math.radians(i)))) + mX
		y2 = (12*chunk*(math.cos(math.radians(i)))) + mY
		## these are broken
		#x2 = (12*ppd*7*(math.cos(math.radians(i)))) + mX
		#y2 = mY - (12*7*ppd*(math.sin(math.radians(i))))  ## because Y is measured from the top left
		line = LineString([(x1,y1),(x2,y2)])
		#print("LINE=",line)
		#print("POLY=",polygon)
		inter = polygon.intersection(line)
		if (isinstance(inter,LineString)):
			#print("INTER=", inter)
			#print("COORDS=", inter.coords[1][0])
			isoIntercepts[i][0] = inter.coords[1][0];
			isoIntercepts[i][1] = inter.coords[1][1];
			#f2.write(str(i) + "\t" + str(inter.coords[1][0]) + "\t" + str(inter.coords[1][1]) + "\n")
		else:
			inter = inter[-1]
			# TODD inter = inter.geoms[-1]
			isoIntercepts[i][0] = inter.coords[1][0];
			isoIntercepts[i][1] = inter.coords[1][1];
			#f2.write(str(i) + "\t" + str(inter.coords[1][0]) + "\t" + str(inter.coords[1][1]) + "\n")

	#print("IIT=",isoIntercepts[1][0])

	to_add = readIntersectFile(mX, mY, ppd, isoIntercepts, new_normData)
	# TODD print (to_add)
	## append to the annotation file
	#fh = open(annot_file, "a")
	f2.write(to_add)

	f2.write(remainder)
	#f2.write("\t</annotations>\n")
	#f2.write("</annotation-per-image>\n")
	#f2.write("</set>\n")
	f2.close()


def usage():
	print("Usage: %s [-i <FILE>]" % os.path.basename(__file__))


if __name__ == "__main__":
	main(sys.argv[1:])
