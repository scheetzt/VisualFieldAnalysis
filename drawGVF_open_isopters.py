#!/usr/bin/python

## Converts Truthmarker generated XML files into pdf figures
## Usage: python3 drawGVF.py [options] [source]
## Options:
##   -i ..., --input=...             Use specified XML file.
##   -h                              print out a help message describing the options
##   -o DIR, --output-dir DIR        specify the directory in which output will be rendered
##   -output-full-path PATH          specifies the filename of the PDF to be rendered. (requires XML has a single field)
##   --no-log                        Do not print out logging messages. Useful when running as a server process.

## DOC STRING attempt
"""Usage: drawGVF.py [-h] [-i FILE | --input FILE] [-o DIR | --output-dir DIR] [--output-file-path FULLPATH] [--no-log] """


from pyx import *
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely import geometry
from shapely.geometry import LineString
from shapely.geometry import LinearRing
from shapely.validation import make_valid
import math
import os, sys, getopt, pprint, re, math
import xml.dom.minidom
import numpy as np
import os
import datetime
import string
from pprint import pprint

def get_imagepath(each, output_dir, output_filepath):
	if(output_filepath != ''):
		return output_filepath
	imagename = each.getElementsByTagName("imagename")
	fileName = imagename[0].firstChild.nodeValue
	fileName = fileName.split('/')
	if len(fileName) > 1:
		fileName = fileName[-1]
	else:	
		fileName = fileName[0]
	return(output_dir + "/" + fileName)

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
	if (distance < 0):
		size = distance/300
		scale = 25*size*-1

	return(scale,num_fid,middleX,middleY,topX,topY,rightX,rightY,bottomX,bottomY,leftX,leftY,eye)

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
		fidComplete.append(xTest)
		fidComplete.append(yTest)
		ptsComplete.append(rFidX_pt[i])
		ptsComplete.append(rFidY_pt[i])
	return (fidComplete, ptsComplete)

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


def collectIsopters(each,avgRotation,avgRRatio,mX,mY):
	annotations = each.getElementsByTagName("annotation")
	isopters = []
	scotomas = []
	openIso = []
	exIso = []
	
	numI1e = 1
	numI2e = 1
	numI3e = 1
	numI4e = 1
	numIII4e = 1
	numV4e = 1

	numI1eSC = 1
	numI2eSC = 1
	numI3eSC = 1
	numI4eSC = 1
	numIII4eSC = 1
	numV4eSC = 1

	numI1eOpen = 1
	numI2eOpen = 1
	numI3eOpen = 1
	numI4eOpen = 1
	numIII4eOpen = 1
	numV4eOpen = 1
	
	numV4eEx = 1

	#Isopters
	for annotation in annotations:
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I1e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				isopters.append('I1e-' + str(numI1e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				numI1e = numI1e + 1
			if (name == 'I1e-sc'):				
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				scotomas.append('I1e-sc' + str(numI1eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				numI1eSC = numI1eSC + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I1e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				openIso.append('I1e-open'+str(numI1eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				numI1eOpen = numI1eOpen + 1
	for annotation in annotations:
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name =='I2e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				isopters.append('I2e-'+ str(numI2e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				numI2e = numI2e + 1
			if (name == 'I2e-sc'):				
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				scotomas.append('I2e-sc' + str(numI2eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				numI2eSC = numI2eSC + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name =='I2e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				openIso.append('I2e-open' + str(numI2eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				numI2eOpen = numI2eOpen + 1
	for annotation in annotations:
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I3e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				isopters.append('I3e-'+str(numI3e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				numI3e = numI3e + 1
			if (name == 'I3e-sc'):				
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				scotomas.append('I3e-sc' + str(numI3eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				numI3eSC = numI3eSC + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I3e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				openIso.append('I3e-open' + str(numI3eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				numI3eOpen = numI3eOpen + 1
	for annotation in annotations:
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I4e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				isopters.append('I4e-' + str(numI4e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				numI4e = numI4e + 1
			if (name == 'I4e-sc'):				
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				scotomas.append('I4e-sc' + str(numI4eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				numI4eSC = numI4eSC + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I4e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				openIso.append('I4e-open' + str(numI4eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				numI4eOpen = numI4eOpen + 1
	for annotation in annotations:
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'III4e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				isopters.append('III4e-' + str(numIII4e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				numIII4e = numIII4e + 1
			if (name == 'III4e-sc'):				
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				scotomas.append('III4e-sc' + str(numIII4eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				numIII4eSC = numIII4eSC + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'III4e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				openIso.append('III4e-open' + str(numIII4eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				numIII4eOpen = numIII4eOpen + 1
	for annotation in annotations:
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'V4e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				isopters.append('V4e-' + str(numV4e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				numV4e = numV4e + 1
			if (name == 'V4e-sc'):				
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				scotomas.append('V4e-sc' + str(numV4eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				numV4eSC = numV4eSC + 1
			if (name == 'V4e-Ex'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				exIso.append('V4e-Ex' + str(numV4eEx))
				exIso.append(xlist2)
				exIso.append(ylist2)
				numV4eEx = numV4eEx + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name =='V4e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				openIso.append('V4e-open' + str(numV4eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				numV4eOpen = numV4eOpen + 1

	return(isopters, scotomas, openIso, exIso)


def construct_hos(isopters, scotomas, openIso, mX, mY, tX, tY, bX, bY, rX, rY, lX, lY, degree_len):

	alphabet_string = string.ascii_lowercase
	polygons_iso = []
	polygons_sc = []
	polygons_open = []

	#isopters and scotomas to Polygons
	for i in range(0,len(isopters),3):
		point_list = []
		for j in range(0,len(isopters[i+1])):
			point_list.append((isopters[i+1][j], isopters[i+2][j]))
		polygon = Polygon(point_list)
		polygons_iso.append(isopters[i])
		polygons_iso.append(polygon)
	for i in range(0,len(scotomas),3):
		point_list = []
		for j in range(0,len(scotomas[i+1])):
			point_list.append((scotomas[i+1][j], scotomas[i+2][j]))
		scotoma = Polygon(point_list)
		polygons_sc.append(scotomas[i])
		polygons_sc.append(scotoma)
	for i in range(0,len(openIso),3):
		point_list = []
		for j in range(0,len(openIso[i+1])):
			point_list.append((openIso[i+1][j], openIso[i+2][j]))
		open_line = LineString(point_list)
		polygons_open.append(openIso[i])
		polygons_open.append(open_line)
		
	all_isopters = polygons_iso + polygons_sc + polygons_open

	isopter_labels_dict = {}
	letters_count = 0
	for i in range(0,len(all_isopters),2):
		letter_name = alphabet_string[letters_count]
		letters_count = letters_count + 1
		isopter_labels_dict[letter_name] = {}
		isopter_labels_dict[letter_name]['name'] = all_isopters[i]
		isopter_labels_dict[letter_name]['polygon'] = all_isopters[i+1]
		if "I1e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
			if "open" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "open"
			isopter_labels_dict[letter_name]['color'] = color.rgb.green
		if "I2e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
			if "open" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "open"
			isopter_labels_dict[letter_name]['color'] = color.rgb.red
		if "I3e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
			if "open" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "open"
			isopter_labels_dict[letter_name]['color'] = color.rgb(.4,.4,.4)
		if "I4e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
			if "open" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "open"
			isopter_labels_dict[letter_name]['color'] = color.rgb.blue
		if "III4e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
			if "open" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "open"
			isopter_labels_dict[letter_name]['color'] = color.rgb(.8,.4,.1)
		if "V4e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
			if "open" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "open"
			isopter_labels_dict[letter_name]['color'] =  color.rgb(1,.3,1)


	
	for each in isopter_labels_dict:
		points = []
		if isopter_labels_dict[each]['type'] == 'open':
			points = isopter_labels_dict[each]['polygon'].coords[:]
		else:
			points = isopter_labels_dict[each]['polygon'].exterior.coords[:]

		ring = LinearRing(points)
		if (ring.is_ccw):
			new_points = points[::-1]
			if (isopter_labels_dict[each]['type'] == 'open'):
				isopter_labels_dict[each]['polygon'] = LineString(new_points)
			else:
				isopter_labels_dict[each]['polygon'] = Polygon(new_points)

	#Find all intersections
	intersects_dict = {}
	for each in isopter_labels_dict:
		intersect_list_sc = []
		intersect_list_iso = []
		intersects_dict[each] = {}
		intersects_dict[each]["isopters"] = {}
		intersects_dict[each]["scotomas"] = {}
		for other in isopter_labels_dict:
			if (each != other and isopter_labels_dict[each]['type'] != 'open' and isopter_labels_dict[other]['type'] != 'open'):
				inter = isopter_labels_dict[each]['polygon'].intersects(isopter_labels_dict[other]['polygon'])
				if (inter == True):
					#find which isopter is inside the other
					i_point = isopter_labels_dict[other]['polygon'].exterior.coords[1]
					i_Point = Point(i_point[0], i_point[1])
					point_test  = isopter_labels_dict[each]['polygon'].contains(i_Point)
					if (point_test == True):
						if ("-sc" in isopter_labels_dict[other]['name']):
							intersect_list_sc.append(other)
						else:
							intersect_list_iso.append(other)
		intersects_dict[each]["isopters"] = intersect_list_iso
		intersects_dict[each]["scotomas"] = intersect_list_sc

	#Prune HOS
	for each in intersects_dict:
		all_children = intersects_dict[each]["isopters"] + intersects_dict[each]["scotomas"]
		sc_to_del = []
		iso_to_del = []
		for each_child in all_children:
			for c2_iso in intersects_dict[each_child]["isopters"]:
				if (c2_iso in intersects_dict[each]["isopters"]):
					iso_to_del.append(c2_iso)
			for c2_sc in intersects_dict[each_child]["scotomas"]:
				if (c2_sc in intersects_dict[each]["scotomas"]):
					sc_to_del.append(c2_sc)
		for iso in iso_to_del:
			if (iso in intersects_dict[each]["isopters"]):
				intersects_dict[each]["isopters"].remove(iso)
		for sc in sc_to_del:
			if (sc in intersects_dict[each]["scotomas"]):
				intersects_dict[each]["scotomas"].remove(sc)

	

	open_isopter_dict = {}
	for each in isopter_labels_dict:
		if (isopter_labels_dict[each]['type'] == 'open'):
			point_test = ''
			bounding_iso = ''
			open_isopter_dict[each] = {}
			parents_list = []
			for other in isopter_labels_dict:
				if (each != other):
					line_point = isopter_labels_dict[each]['polygon'].coords[0]
					line_Point = Point(line_point[0], line_point[1])
					point_test  = isopter_labels_dict[other]['polygon'].contains(line_Point)
					if (point_test == True):	
						parents_list.append(other)
			open_isopter_dict[each] = parents_list

	#Make root list
	root_list = []
	for each in isopter_labels_dict:
		root_list.append(each)

	for each in isopter_labels_dict:
		for other in intersects_dict:
			if (each in intersects_dict[other]['isopters']):
				if (each in root_list):
					root_list.remove(each)
			if (each in intersects_dict[other]['scotomas']):
				if (each in root_list):
					root_list.remove(each)

	#Works for many open isopters as long as the isopters are "in order" from smallest to largest
	for each in open_isopter_dict:
		if (len(open_isopter_dict[each]) == 0):
			intersects_dict[each]["isopters"] = [root_list[0]]
			intersects_dict[each]["scotomas"] = []
			root_list = [each]
		#elif (len(open_isopter_dict[each]) == 1):
		else:
			new_parent = open_isopter_dict[each][0]
			old_child_iso = intersects_dict[new_parent]['isopters']
			old_child_sc = intersects_dict[new_parent]['scotomas']
			iso_list = []
			iso_list.append(each)
			intersects_dict[new_parent]['isopters'] = iso_list
			intersects_dict[each]['isopters'] = old_child_iso
			intersects_dict[each]['scotomas'] = old_child_sc

	#Need to make object that contains the parent of each entry
	parent_dict = {}
	for each in isopter_labels_dict:
		parent_dict[each] = {}
		for other in intersects_dict:
			if each in intersects_dict[other]["isopters"]:
				parent_dict[each] = other
			if each in intersects_dict[other]["scotomas"]:
				parent_dict[each] = other

	new_points = []
	double_open_iso = 0
	for each in open_isopter_dict:
		iso_name = isopter_labels_dict[each]['name']
		if '-open2' in iso_name:
			double_open_iso = 1
		new_points.append(each)
		new_x, new_y = complete_isopter(each, isopter_labels_dict, intersects_dict, root_list, parent_dict, new_points, degree_len)
		new_points.append(new_x)
		new_points.append(new_y)

	#print ("new_points")
	#print (new_points)

	return (isopter_labels_dict, intersects_dict, root_list, parent_dict, new_points, double_open_iso)

def test_orth(isopter_labels_dict, child_seg_1, child_seg_2, iso_extrap, degree_len):

	j = 0
	empty_count_seg_1 = 0
	point_count_seg_1 = 0
	empty_count_seg_2 = 0
	point_count_seg_2 = 0

	for i in range(1,(len(child_seg_1)-1),3):
		x1 = child_seg_1[i-1][0]
		y1 = child_seg_1[i-1][1]
		x2 = child_seg_1[i+1][0]
		y2 = child_seg_1[i+1][1]
		x3 = child_seg_1[i][0]
		y3 = child_seg_1[i][1]
		delta_x = (x2-x1)
		delta_y = (y2-y1)
		if (delta_y == 0):
			delta_y = 0.1
		orth_slope = (delta_x/delta_y)*(-1)
		r = math.sqrt(1+(orth_slope*orth_slope))
		proj_x = 0
		proj_y = 0
		total_projection_distance = degree_len*90
		if (delta_y > 0):
			proj_y = y3 - ((total_projection_distance*orth_slope)/r)
			proj_x = x3 - (total_projection_distance/r)
		else:
			proj_y = y3 + ((total_projection_distance*orth_slope)/r)
			proj_x = x3 + (total_projection_distance/r)
		projection_line = LineString([(x3,y3), (proj_x, proj_y)])
		intersect_point_open = projection_line.intersection(isopter_labels_dict[iso_extrap]['polygon'])
		if (intersect_point_open.geom_type == "Point"):
			point_count_seg_1 = point_count_seg_1 + 1
		elif (intersect_point_open.geom_type == "MultiPoint"):
			point_count_seg_1 = point_count_seg_1 + 1
		else:
			empty_count_seg_1 = empty_count_seg_1 + 1
		j = j + 3

	for i in range(1,(len(child_seg_2)-1),3):
		x1 = child_seg_2[i-1][0]
		y1 = child_seg_2[i-1][1]
		x2 = child_seg_2[i+1][0]
		y2 = child_seg_2[i+1][1]
		x3 = child_seg_2[i][0]
		y3 = child_seg_2[i][1]
		delta_x = (x2-x1)
		delta_y = (y2-y1)
		if (delta_y == 0):
			delta_y = 0.1
		orth_slope = (delta_x/delta_y)*(-1)
		r = math.sqrt(1+(orth_slope*orth_slope))
		proj_x = 0
		proj_y = 0
		total_projection_distance = degree_len*90
		if (delta_y > 0):
			proj_y = y3 - ((total_projection_distance*orth_slope)/r)
			proj_x = x3 - (total_projection_distance/r)
		else:
			proj_y = y3 + ((total_projection_distance*orth_slope)/r)
			proj_x = x3 + (total_projection_distance/r)
		projection_line = LineString([(x3,y3), (proj_x, proj_y)])
		intersect_point_open = projection_line.intersection(isopter_labels_dict[iso_extrap]['polygon'])
		if (intersect_point_open.geom_type == "Point"):
			point_count_seg_2 = point_count_seg_2 + 1
		elif (intersect_point_open.geom_type == "MultiPoint"):
			point_count_seg_2 = point_count_seg_2 + 1
		else:
			empty_count_seg_2 = empty_count_seg_2 + 1
		j = j + 3

	if (point_count_seg_1 == 0):
		point_count_seg_1 = 0.1
	if (point_count_seg_2 == 0):
		point_count_seg_2 = 0.1

	seg_1_ratio = float(empty_count_seg_1/point_count_seg_1)
	seg_2_ratio = float(empty_count_seg_2/point_count_seg_2)

	#going to want to keep the segment with the higher ratio (more empty points)
	if (seg_1_ratio > seg_2_ratio):
		return (child_seg_1)
	else:
		return (child_seg_2)
		

	#print ("seg 1 ratio: ", seg_1_ratio)
	#print ("seg 2 ratio: ", seg_2_ratio)

	#print ("point count seg 1: ", point_count_seg_1)
	#print ("empty count seg 1: ", empty_count_seg_1)
	#print ("point count seg 2: ", point_count_seg_2)
	#print ("empty count seg 2: ", empty_count_seg_2)
	

def complete_isopter(each, isopter_labels_dict, intersects_dict, root_list, parent_dict, new_points, degree_len):

	iso_extrap = each
	iso_child = intersects_dict[each]['isopters'][0]
#	if (isopter_labels_dict[iso_child]['type'] == 'open'):
#		current_child = iso_child
#		while (isopter_labels_dict[current_child]['type'] == 'open'):
#			current_child = intersects_dict[current_child]['isopters'][0]
#			iso_child = current_child
	iso_parent = parent_dict[each]
	if (iso_parent != {}):
		if (isopter_labels_dict[iso_parent]['type'] == 'open'):
			iso_parent = {}

	iso_extrap = each
	print (isopter_labels_dict[each]['name'])

	first_point = isopter_labels_dict[each]['polygon'].coords[0]
	last_point = isopter_labels_dict[each]['polygon'].coords[-1]
	
	points_list_iso = isopter_labels_dict[each]['polygon'].coords

	if (isopter_labels_dict[iso_child]['type'] == 'open'):
		combined_poly_points = []
		points_list_child_1 = isopter_labels_dict[iso_child]['polygon'].coords[:]
		for each in points_list_child_1:
			#print (each)
			combined_poly_points.append(each)
		new_open_line = ''
		#print (points_list_child)
		iso_name  = isopter_labels_dict[iso_child]['name']
		for i in range(0,len(new_points),3):
			test_name = isopter_labels_dict[new_points[i]]['name']
			if (test_name == isopter_labels_dict[iso_child]['name']):
				to_combine_points = []
				for j in range(0,len(new_points[i+1])):
					to_combine_points.append((new_points[i+1][j], new_points[i+2][j]))
				new_open_line = LineString(to_combine_points)
		points_list_child_2 = new_open_line.coords[:]
		for each in points_list_child_2:
			#print (each)
			combined_poly_points.append(each)
		full_poly = Polygon(combined_poly_points)
		points_list_child = full_poly.exterior.coords
			
	else:
		points_list_child = isopter_labels_dict[iso_child]['polygon'].exterior.coords

	#Need to find closest point on full polygon
	# For the child isopter
	cfp_1, fp_index_1, min_fp_dist_1 = closest_to_first_point(points_list_iso,points_list_child, first_point)
	clp_1, lp_index_1, min_lp_dist_1 = closest_to_last_point(points_list_iso,points_list_child, last_point)
	
	if (iso_parent != {}):
		points_list_parent = isopter_labels_dict[iso_parent]['polygon'].exterior.coords
		cfp_2, fp_index_2, min_fp_dist_2 = closest_to_first_point(points_list_iso,points_list_parent, first_point)
		clp_2, lp_index_2, min_lp_dist_2 = closest_to_last_point(points_list_iso,points_list_parent, last_point)
	else: 
		cfp_2 = 0
		fp_index_2 = 0
		min_fp_dist_2 = 0
		clp_2 = 0
		lp_index_2 = 0
		min_lp_dist_2 = 0

	if (fp_index_1 < lp_index_1):
		print ("case 1")
		child_seg_1 = points_list_child[(fp_index_1):(lp_index_1)]
		child_seg_2 = points_list_child[(lp_index_1):] + points_list_child[:(fp_index_1)]

		#print ("segment 1")
		#for each in child_seg_1:
		#	print (each)
		#print ("segment 2")
		#for each in child_seg_2:
		#	print (each)

#####Going to add a function to test the orthgonal line off of the segments to see if they cross the open isopter
		seg_to_extrap = test_orth(isopter_labels_dict, child_seg_1, child_seg_2, iso_extrap, degree_len)
		#for each in seg_to_extrap:
		#	print (each)
		#print (seg_to_extrap)
		new_x, new_y = extrap_points(seg_to_extrap, iso_parent, isopter_labels_dict, min_fp_dist_1, min_fp_dist_2, min_lp_dist_1, min_lp_dist_2)
		print ("extrapolated: ")
		for i in range(0,len(new_x)):
			print (new_x[i], new_y[i])
		return (new_x, new_y)

#		to_add_x_1 = []
#		to_add_y_1 = []
#		#to_add_x_1.append(last_point[0])
#		#to_add_y_1.append(last_point[1]))
#		new_x, new_y = extrap_points(child_seg_1, iso_parent, isopter_labels_dict, min_fp_dist_1, min_fp_dist_2, min_lp_dist_1, min_lp_dist_2)
#		#print ("Extrapolated points")
#		#print ("x")
#		for each in new_x:
#			to_add_x_1.append(each)
#		#	print (each)
#		#print ("y")
#		for each in new_y:
#			to_add_y_1.append(each)
#		#	print (each)
#		#to_add_x_1.append(first_point[0])
#		#to_add_y_1.append(first_point[1])
#
#		to_add_x_2 = []
#		to_add_y_2 = []
#		#to_add_x_2.append(last_point[0])
#		#to_add_y_2.append(last_point[1])
#		new_x, new_y = extrap_points(child_seg_2, iso_parent, isopter_labels_dict, min_fp_dist_1, min_fp_dist_2, min_lp_dist_1, min_lp_dist_2)
#		#print ("Extrapolated points:")
#		#print ("x")
#		for each in new_x:
#			to_add_x_2.append(each)
#		#	print (each)
#		#print ("y")
#		for each in new_y:
#			to_add_y_2.append(each)
#		#	print (each)
#		#to_add_x_2.append(first_point[0])
#		#to_add_y_2.append(first_point[1])
#
#		#Find which extrapolated segment is closer to the open isopter
#
#		#Average distance to segment 1
#		total_distances_1 = 0
#		total_points_1 = len(to_add_x_1)
#		average_distance_1 = 0
#
#		for i in range (0, len(to_add_x_1)):
#			min_distance = 10000
#			for open_point in points_list_iso:
#				current_distance = math.hypot(open_point[0]-to_add_x_1[i],open_point[1]-to_add_y_1[i])
#				if (current_distance < min_distance):
#					min_distance = current_distance
#			total_distances_1 = total_distances_1 + min_distance
#		average_distance_1 = total_distances_1/total_points_1
#		
#		#Average distance to segment 2
#		total_distances_2 = 0
#		total_points_2 = len(to_add_x_2)
#		average_distance_2 = 0
#
#		for i in range (0,len(to_add_x_2)):
#			min_distance = 10000
#			for open_point in points_list_iso:
#				current_distance = math.hypot(open_point[0]-to_add_x_2[i],open_point[1]-to_add_y_2[i])
#				if (current_distance < min_distance):
#					min_distance = current_distance
#			total_distances_2 = total_distances_2 + min_distance
#
#		average_distance_2 = total_distances_2/total_points_2
#
#		#if (average_distance_2 > average_distance_1):
#		#	check_linearRing(to_add_x_2, to_add_y_2)
#		#	return (to_add_x_2, to_add_y_2)
#		#else:
#		#	check_linearRing(to_add_x_1, to_add_y_1)
#		#	return (to_add_x_1, to_add_y_1)
#
#		print ("From linear ring:")
#		if (average_distance_2 > average_distance_1):
#			print ("returning segment 2")
#			to_add_x, to_add_y = check_linearRing(to_add_x_2, to_add_y_2)
#			for i in range (0,len(to_add_x)):
#				print (to_add_x[i], to_add_y[i])
#			return (to_add_x, to_add_y)
#		else:
#			print ("returning segment 1")
#			to_add_x, to_add_y = check_linearRing(to_add_x_1, to_add_y_1)
#			#for i in range (0,len(to_add_x)):
#			#	print (to_add_x[i], to_add_y[i])
#			return (to_add_x, to_add_y)
	else:
		print ("case 2")
		child_seg_1 = points_list_child[(fp_index_1):] + points_list_child[:lp_index_1]
		child_seg_2 = points_list_child[(lp_index_1):(fp_index_1)]

		#print ("segment 1")
		#for each in child_seg_1:
		#	print (each)
		#print ("segment 2")
		#for each in child_seg_2:
		#	print (each)

		seg_to_extrap = test_orth(isopter_labels_dict, child_seg_1, child_seg_2, iso_extrap, degree_len)
		new_x, new_y = extrap_points(seg_to_extrap, iso_parent, isopter_labels_dict, min_fp_dist_1, min_fp_dist_2, min_lp_dist_1, min_lp_dist_2)
		return (new_x, new_y)

#		to_add_x_1 = []
#		to_add_y_1 = []
#		#to_add_x_1.append(first_point[0])
#		#to_add_y_1.append(first_point[1])
#		new_x, new_y = extrap_points(child_seg_1, iso_parent, isopter_labels_dict, min_fp_dist_1, min_fp_dist_2, min_lp_dist_1, min_lp_dist_2)
#		for each in new_x:
#			to_add_x_1.append(each)
#		for each in new_y:
#			to_add_y_1.append(each)
#		#to_add_x_1.append(last_point[0])
#		#to_add_y_1.append(last_point[1])
#	
#		to_add_x_2 = []
#		to_add_y_2 = []
#		#to_add_x_2.append(first_point[0])
#		#to_add_y_2.append(first_point[1])
#		new_x, new_y = extrap_points(child_seg_2, iso_parent, isopter_labels_dict, min_fp_dist_1, min_fp_dist_2, min_lp_dist_1, min_lp_dist_2)
#		for each in new_x:
#			to_add_x_2.append(each)
#		for each in new_y:
#			to_add_y_2.append(each)
#		#to_add_x_2.append(last_point[0])
#		#to_add_y_2.append(last_point[1])
#
#		#Find which extrapolated segment is closer to the open isopter
#
#		#Average distance to segment 1
#		total_distances_1 = 0
#		total_points_1 = len(to_add_x_1)
#		average_distance_1 = 0
#
#		for i in range (0, len(to_add_x_1)):
#			min_distance = 10000
#			for open_point in points_list_iso:
#				current_distance = math.hypot(open_point[0]-to_add_x_1[i],open_point[1]-to_add_y_1[i])
#				if (current_distance < min_distance):
#					min_distance = current_distance
#			total_distances_1 = total_distances_1 + min_distance
#		average_distance_1 = total_distances_1/total_points_1
#		
#		#Average distance to segment 2
#		total_distances_2 = 0
#		total_points_2 = len(to_add_x_2)
#		average_distance_2 = 0
#
#		for i in range (0,len(to_add_x_2)):
#			min_distance = 10000
#			for open_point in points_list_iso:
#				current_distance = math.hypot(open_point[0]-to_add_x_2[i],open_point[1]-to_add_y_2[i])
#				if (current_distance < min_distance):
#					min_distance = current_distance
#			total_distances_2 = total_distances_2 + min_distance
#
#		average_distance_2 = total_distances_2/total_points_2
#
#		#print ("average_distance_1, average_distance_2")
#		#print (average_distance_1, average_distance_2)

#		if (average_distance_2 > average_distance_1):
#			to_add_x, to_add_y = check_linearRing(to_add_x_2, to_add_y_2)
#			#for i in range (0,len(to_add_x)):
#			#	print (to_add_x[i], to_add_y[i])
#			return (to_add_x, to_add_y)
#		else:
#			to_add_x, to_add_y = check_linearRing(to_add_x_1, to_add_y_1)
#			#for i in range (0,len(to_add_x)):
#			#	print (to_add_x[i], to_add_y[i])
#			return (to_add_x, to_add_y)

def check_linearRing(to_add_x, to_add_y):
	points = []
	#print ("input points")
	for i in range(0, len(to_add_x)):
		points.append([to_add_x[i], to_add_y[i]])
	#	print(to_add_x[i], to_add_y[i])
	original_ring = Polygon(points)
	#print ("checking")

	if (original_ring.is_valid == False):
		all_poly = make_valid(original_ring)
		largest_len  = 0
		largest_poly = 0
		#print ("invalid polygon")
		new_x = []
		new_y = []

		if (all_poly.geom_type == 'MultiPolygon'):	
			#print ("I'm a multipolygon")
			for each in all_poly:
				#print ("poly_chunk")
				#print (each)
				#print (each.length)
				#print (each.exterior.coords[:])
				#print (len(each.geoms))
				current_points = each.exterior.coords[:]
				#for point in current_points:
				#	print (point)
				if (len(current_points) > largest_len):
					largest_len = len(current_points)
					largest_poly = each
			#print ("largest_poly")
			#print (largest_poly)
			#print ("ext coords")
			#print (largest_poly.exterior.coords[:])
			#Kind of a stab in the dark here with the trimming 
			new_points = largest_poly.exterior.coords[1:-4]
			for each in new_points:
				new_x.append(each[0])
				new_y.append(each[1])


		elif (all_poly.geom_type == 'Polygon'):	
			current_points = all_poly.exterior.coords[:-2]
			for each in current_points:
				new_x.append(each[0])
				new_y.append(each[1])

#		else:
#			#print ("Other")


		#return (new_x, new_y)
		return (to_add_x, to_add_y)

	else:
		#print ("is valid")
		return (to_add_x, to_add_y)
		

	
def extrap_points(points_to_extrap, iso_parent, isopter_labels_dict, min_fp_dist_1, min_fp_dist_2, min_lp_dist_1, min_lp_dist_2):


	to_add_x = []
	to_add_y = []
	if (iso_parent == {}):
		distance_difference = min_fp_dist_1 - min_lp_dist_1
		distance_increment = (distance_difference/len(points_to_extrap))
		#print (distance_difference, distance_increment)
		j = 0
		for i in range(1,(len(points_to_extrap)-1)):
			x1 = points_to_extrap[i-1][0]
			y1 = points_to_extrap[i-1][1]
			x2 = points_to_extrap[i+1][0]
			y2 = points_to_extrap[i+1][1]
			x3 = points_to_extrap[i][0]
			y3 = points_to_extrap[i][1]
			delta_x = (x2-x1)
			delta_y = (y2-y1)
			if (delta_y == 0):
				delta_y = 0.1
			orth_slope = (delta_x/delta_y)*(-1)
			r = math.sqrt(1+(orth_slope*orth_slope))
			distance = min_lp_dist_1+(distance_increment*j)
			new_x = 0
			new_y = 0
			if (delta_y > 0):
				new_y = y3 - ((distance*orth_slope)/r)
				new_x = x3 - (distance/r)
			else:
				new_y = y3 + ((distance*orth_slope)/r)
				new_x = x3 + (distance/r)
			to_add_x.append(new_x)
			to_add_y.append(new_y)
			j = j+1
	else:
		rel_start = (min_lp_dist_1/(min_lp_dist_1 + min_lp_dist_2))
		rel_end = (min_fp_dist_1/(min_fp_dist_1 + min_fp_dist_2))
		rel_increment = (rel_end - rel_start) / (len(points_to_extrap))
		fp_distance = min_fp_dist_1+min_fp_dist_2
		lp_distance = min_lp_dist_1+min_lp_dist_2
		total_projection_distance = lp_distance * 5

		j = 0
		for i in range(1,(len(points_to_extrap)-1),4):
			x1 = points_to_extrap[i-1][0]
			y1 = points_to_extrap[i-1][1]
			x2 = points_to_extrap[i+1][0]
			y2 = points_to_extrap[i+1][1]
			x3 = points_to_extrap[i][0]
			y3 = points_to_extrap[i][1]
			delta_x = (x2-x1)
			delta_y = (y2-y1)
			if (delta_y == 0):
				delta_y = 0.1
			orth_slope = (delta_x/delta_y)*(-1)
			r = math.sqrt(1+(orth_slope*orth_slope))

			proj_x = 0
			proj_y = 0
			if (delta_y > 0):
				proj_y = y3 - ((total_projection_distance*orth_slope)/r)
				proj_x = x3 - (total_projection_distance/r)
			else:
				proj_y = y3 + ((total_projection_distance*orth_slope)/r)
				proj_x = x3 + (total_projection_distance/r)

			projection_line = LineString([(x3,y3), (proj_x, proj_y)])
			intersect_point = projection_line.intersection(isopter_labels_dict[iso_parent]['polygon'])
			if (intersect_point.geom_type == "MultiLineString"):
				intersect_point = intersect_point.geoms[0]
			projection_distance = math.hypot(intersect_point.coords[0][0]-intersect_point.coords[1][0],intersect_point.coords[0][1]-intersect_point.coords[1][1])
			rel_spot = (rel_start + (rel_increment*j))
			rel_dist = rel_spot * projection_distance

			new_x = 0
			new_y = 0
			if (delta_y > 0):
				new_y = y3 - ((rel_dist*orth_slope)/r)
				new_x = x3 - (rel_dist/r)
			else:
				new_y = y3 + ((rel_dist*orth_slope)/r)
				new_x = x3 + (rel_dist/r)
			to_add_x.append(new_x)
			to_add_y.append(new_y)
			#print (new_x, new_y)
			#print (each)
			j = j+4

	return(to_add_x,to_add_y)

def closest_to_first_point(points_list_iso,points_list, first_point):
	
	min_first_point_dist = 10000
	closest_first_point = []
	before_closest_first_point = []
	after_closest_first_point = []
	first_point_index = 0 

	#Find point closest to first point
	for i in range(0,len(points_list)):
		current_distance = math.hypot(points_list[i][0]-first_point[0],points_list[i][1]-first_point[1])
		if (current_distance < min_first_point_dist):
			min_first_point_dist = current_distance
			before_closest_first_point = [points_list[i-1][0],points_list[i-1][1]]
			closest_first_point = [points_list[i][0],points_list[i][1]]
			first_point_index = i + 1
			after_closest_first_point = [points_list[i+1][0],points_list[i+1][1]]
	before_first_min_dist_open = 10000
	after_first_min_dist_open = 1000
	for i in range(0, len(points_list_iso)):
		current_dist_before = math.hypot(points_list_iso[i][0]-before_closest_first_point[0],points_list_iso[i][1]-before_closest_first_point[1])
		current_dist_after = math.hypot(points_list_iso[i][0]-after_closest_first_point[0],points_list_iso[i][1]-after_closest_first_point[1])
		if (current_dist_before < before_first_min_dist_open):
			before_first_min_dist_open = current_dist_before
		if (current_dist_after < after_first_min_dist_open):
			after_first_min_dist_open = current_dist_after
	return (closest_first_point, first_point_index, min_first_point_dist)

def closest_to_last_point(points_list_iso, points_list, last_point):
	
	
	min_last_point_dist = 10000
	closest_last_point = []
	before_closest_last_point = []
	after_closest_last_point = []
	last_point_index = 0

	#Find point closest to last point
	for i in range(0,len(points_list)):
		current_distance = math.hypot(points_list[i][0]-last_point[0],points_list[i][1]-last_point[1])
		if (current_distance < min_last_point_dist):
			min_last_point_dist = current_distance
			before_closest_last_point = [points_list[i-1][0],points_list[i-1][1]]
			closest_last_point = [points_list[i][0],points_list[i][1]]
			last_point_index = i - 1
			after_closest_last_point = [points_list[i+1][0],points_list[i+1][1]]
	before_last_min_dist_open = 10000
	after_last_min_dist_open = 1000
	for i in range(0, len(points_list_iso)):
		current_dist_before = math.hypot(points_list_iso[i][0]-before_closest_last_point[0],points_list_iso[i][1]-before_closest_last_point[1])
		current_dist_after = math.hypot(points_list_iso[i][0]-after_closest_last_point[0],points_list_iso[i][1]-after_closest_last_point[1])
		if (current_dist_before < before_last_min_dist_open):
			before_last_min_dist_open = current_dist_before
		if (current_dist_after < after_last_min_dist_open):
			after_last_min_dist_open = current_dist_after
	if (last_point_index < 0):
		last_point_index = 0
	return (closest_last_point, last_point_index, min_last_point_dist)
	
def un_flip_y(y, mY):
	new_y = (y - mY)
	new_y = (new_y)*(-1)
	new_y = (new_y+mY)
	
	return (new_y)
	

def update_xml(new_points, inputfile, isopter_labels_dict, add_to_trace, mY):

	print (len(new_points))
	print (new_points[0])
	print (len(new_points[1]))
	print (new_points[1])
	print (len(new_points[2]))
	print (new_points[2])
	

	if (add_to_trace == 1):
		print (inputfile)
		#print (new_points)
		file_name_arr = inputfile.split(".")
		new_file_name = file_name_arr[0] + ".extrapolated.txt"
		#print (new_file_name)
		f_new = open(new_file_name, "w")
		#with open(inputfile) as file_in:
		#	for line in file_in:
		#		if not ('</annotations>' in line):
		#			f_new.write(line)
		#		else:
		#			break


		remainder = ""
		remainder_flag = 0
		## read-through and re-write all of the lines prior to the end of the annotations
		with open(inputfile) as file_in:
			for line in file_in:
				if remainder_flag == 0:
					if not ('</annotations>' in line):
						line.rstrip("\n")
						# TODD print (line)
						f_new.write(line)
					else:
						## read in the remainder of the file into xml_remainder, to print at the end
						remainder_flag = 1
						remainder = line
				else:
					remainder += line

		for i in range(0, len(new_points)-1, 3):
			iso_letter = new_points[i]
			full_name = isopter_labels_dict[iso_letter]['name']
			name = full_name.split("-")[0]
			new_name = name + "-extrap"
			f_new.write("\t\t<annotation>\n")
			f_new.write("\t\t\t<name>"+new_name+"</name>\n")
			f_new.write("\t\t\t<controltype>SplineRoi</controltype>\n")

			f_new.write("\t\t\t<polygon>\n")
			polygon_points = 0
			for x in range (0,len(new_points[i+1]),5):
				polygon_points = polygon_points + 1
			#polygon_points = int((len(new_points[i+1]))/5)+1
			f_new.write("\t\t\t\t<npoints>"+str(polygon_points)+"</npoints>\n")
			f_new.write("\t\t\t\t<xpoints>\n")
			for x in range (0,len(new_points[i+1]),5):
				f_new.write("\t\t\t\t\t<int>"+str(int(new_points[i+1][x]))+"</int>\n")
			f_new.write("\t\t\t\t</xpoints>\n")
			f_new.write("\t\t\t\t<ypoints>\n")
			for y in range (0,len(new_points[i+2]),5):
				new_y = un_flip_y(new_points[i+2][y], mY)
				f_new.write("\t\t\t\t\t<int>"+str(int(new_y))+"</int>\n")
			f_new.write("\t\t\t\t</ypoints>\n")
			f_new.write("\t\t\t</polygon>\n")

			f_new.write("\t\t\t<spline>\n")
			f_new.write("\t\t\t\t<npoints>"+str(len(new_points[i+1]))+"</npoints>\n")
			f_new.write("\t\t\t\t<xpoints>\n")
			for x in new_points[i+1]:
				f_new.write("\t\t\t\t\t<int>"+str(int(x))+"</int>\n")
			f_new.write("\t\t\t\t</xpoints>\n")
			f_new.write("\t\t\t\t<ypoints>\n")
			for y in new_points[i+2]:
				y = un_flip_y(y, mY)
				f_new.write("\t\t\t\t\t<int>"+str(int(y))+"</int>\n")
			f_new.write("\t\t\t\t</ypoints>\n")
			f_new.write("\t\t\t</spline>\n")
			f_new.write("\t\t</annotation>\n")
		f_new.write(remainder)
		#f_new.write("\t</annotations>\n")
		#f_new.write("</annotation-per-image>\n")
		## TODD - This isn't needed in the ABCA4 dataset, and is actively breaking things
		#f_new.write("</set>\n")
		f_new.close()
	else:
		print (inputfile)
		#print (new_points)
		file_name_arr = inputfile.split(".")
		new_file_name = file_name_arr[0] + ".extrapolated.txt"
		#print (new_file_name)
		f_new = open(new_file_name, "w")
		with open(inputfile) as file_in:
			for line in file_in:
				#if not ('</annotations>' in line):
				f_new.write(line)
				#else:
				#	break
		f_new.close()

def parseAnnotation(inputfile, output_dir, output_filepath, use_log):
	inputfile = inputfile
	dom = xml.dom.minidom.parse(inputfile)
	images = dom.getElementsByTagName("annotation-per-image")
	log_fh = None
	logName = "log.txt"

	if (use_log):
		#log_fh = open(logName,"w")
		#log_fh.write("Image" + "\t" + "Sufficient Fiducials" + "\t" + "Marked Curves" + "\t" + "Eye Annotated" + "\t" + "PDF Generated" + "\n")
		for each in images:
			#Name of the output PDF
			fileName = get_imagepath(each, output_dir, output_filepath)
			#Check fiducials and find scaling factor
			scale, num_fid,mX,mY,tX,tY,rX,rY,bX,bY,lX,lY,eye = fiducials_scale(each)
			if (mX == 0 and mY == 0 and scale == 0):
				mX = 10
				mY = 10
				scale = 1
			if (1):
				sufficientFid = "Yes"
				#Calculate average rotation and RRatio for 
				avgRotation, avgRRatio = rotationCalculation(mX,mY,tX,tY,rX,rY,bX,bY,lX,lY)
				#Rotate fiducial points
				rotatedFids, pts = rotateFiducials(avgRotation,avgRRatio,mX,mY,tX,tY,rX,rY,bX,bY,lX,lY)
				rotated_ty = rotatedFids[1]
				chunk_len = (float((mY-rotated_ty)/7))
				degree_len = chunk_len/10
				#Will rotate points and flip y-axis
				isopters, scotomas, openIso, extrapIso = collectIsopters(each,avgRotation,avgRRatio,mX,mY)
				isopter_labels_dict, intersects_dict, root_list, parent_dict, new_points, double_open_iso = construct_hos(isopters, scotomas, openIso, mX, mY, tX, tY, bX, bY, rX, rY, lX, lY, degree_len)
				add_to_trace = 0
				if (double_open_iso == 0 and len(new_points) != 0):
					add_to_trace = 1
				update_xml(new_points, inputfile, isopter_labels_dict, add_to_trace, mY)


def main(argv):
	inputfile = ''
	output_dir = './'
	output_filepath = ''
	use_log = True
	try:
		opts, args = getopt.getopt(argv,"huo:pi:s",["input=", "output-dir=", "output-file-path=", "no-log"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			usage()
			sys.exit()
		elif opt == '--no-log':
			use_log = False
		elif opt in ("--output-file-path"):
			output_filepath = arg
			#print ('Output filepath is ', output_filepath)
		elif opt in ("-o", "--output-dir"):
			output_dir = arg
			#print ('Output Directory is ', output_dir)
			check_linearRing(to_add_x_2, to_add_y_2)
			check_linearRing(to_add_x_2, to_add_y_2)
		elif opt in ("-i", "--input"):
			inputfile = arg
			#print ('Input file is ', inputfile)
			#print ('\n')
	if(inputfile == ''):
		usage()
		sys.exit()
	parseAnnotation(inputfile, output_dir, output_filepath, use_log)

def usage():
	print (__doc__)
if __name__ == "__main__":
	main(sys.argv[1:])
