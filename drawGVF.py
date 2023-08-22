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
	#print("GI:IMAGENAME=", imagename)
	fileName = imagename[0].firstChild.nodeValue
	fileName = fileName.split('/')
	if len(fileName) > 1:
		fileName = fileName[-1]
	else:	
		fileName = fileName[0]
	#print("GI:FILENAME=", fileName)
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
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "90,70    " or annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "90,70"):
			topX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			topY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "270,70    " or annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "270,70"):
			bottomX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			bottomY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "180,90    " or annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "180,90"):
			leftX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			leftY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "0,90    " or annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "0,90"):
			rightX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			rightY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "0,0    " or annotation.getElementsByTagName("name")[0].firstChild.nodeValue == "0,0"):
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
			#distance = abs(middleY - topY)
			#distance_alt = abs(middleX - topX)
		elif (bottomX != 0 and bottomY != 0):
			distance = math.sqrt(((bottomX-middleX)*(bottomX-middleX))+((bottomY-middleY)*(bottomY-middleY)))
			#distance = abs(bottomY - middleY)
			#distance_alt = abs(bottomX - middleX)
		elif (leftX != 0 and leftY != 0):
			distance = math.sqrt(((leftX-middleX)*(leftX-middleX))+((leftY-middleY)*(leftY-middleY)))
			#distance = abs((7/9)*(middleX - leftX))
			#distance_alt = abs((7/9)*(middleY - leftY))
		elif (rightX != 0 and rightY != 0):
			distance = math.sqrt(((rightX-middleX)*(rightX-middleX))+((rightY-middleY)*(rightY-middleY)))
			#distance = abs((7/9)*(rightX - middleX))
			#distance_alt = abs((7/9)*(rightY - middleY))
		else:
			distance = 0
#	if (distance_alt > distance):
#		distance = distance_alt
	if (distance > 0):
		size = distance/300
		scale = 25*size
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
	#print("DEBUG: Length = ", length)
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
		#temptemp = pr/loc[i+1]
		#print("DEBUG: I=", i, " partial_length=", temptemp)
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
	for i in range(len(rFidX)):
		fidComplete.append(rFidX[i])
		fidComplete.append(rFidY[i])
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

## Function to parse the list of isopters tested from the XML tags
def getIsoptersTested(each):
	isoptersTested = []
	isoptersTestedAnnotated = 0
	#to_check = ["I1e", "I2e", "I3e", "I4e", "III4e", "V4e"]
	#for isopter in to_check:
	#	tester = "tested_" + isopter
	#	temp = each.getElementsByTagName(tester)
	#	if temp:
	#		if(temp[0].firstChild.nodeValue == "1"):
	#			isoptersTested.append(isopter)
	#			print(isopter + " present")
	#		else:
	#			print(isopter + " not present")

	tested = each.getElementsByTagName("isopters_tested")
	if(len(tested) > 0):
		isopters = tested[0].childNodes

        	# Parsing the isopters Isopters
		for isopter in isopters:
			if(isopter.nodeType == 1):
				isoptersTestedAnnotated = 1	## annotations are found
				if(isopter.firstChild.nodeValue == "tested"):
					isoptersTested.append(isopter.nodeName.replace("tested_",""))


	#print(isoptersTested)
	## if not annotated, infer isopters based upon those that have been traced

	isoptersInferred = {}
	inferredArray = []

	annotations = each.getElementsByTagName("annotation")

	#Isopters
	for annotation in annotations:
		name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
		name = name.replace(" ","")
		if (name == 'V4e' or name == 'V4e-sc' or name == 'V4e-open'):
			isoptersInferred['V4e'] = 1
		elif (name == 'III4e' or name == 'III4e-sc' or name == 'III4e-open'):
			isoptersInferred['III4e'] = 1
		elif (name == 'I4e' or name == 'I4e-sc' or name == 'I4e-open'):
			isoptersInferred['I4e'] = 1
		elif (name == 'I3e' or name == 'I3e-sc' or name == 'I3e-open'):
			isoptersInferred['I3e'] = 1
		elif (name == 'I2e' or name == 'I2e-sc' or name == 'I2e-open'):
			isoptersInferred['I2e'] = 1
		elif (name == 'I1e' or name == 'I1e-sc' or name == 'I1e-open'):
			isoptersInferred['I1e'] = 1

	#print(isoptersInferred)
	for key in isoptersInferred:
		inferredArray.append(key)
	#print(inferredArray)
	return(isoptersTested, inferredArray, isoptersTestedAnnotated)

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
	exIsoFull = []
	
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
	
	numI3eEx = 1
	numIII4eEx = 1
	numV4eEx = 1

	numIII4eExFull = 1
	numV4eExFull = 1

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
			if (name == 'I3e-extrap'):			
				print ("I3e")	
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				exIso.append('I3e-Ex' + str(numI3eEx))
				exIso.append(xlist2)
				exIso.append(ylist2)
				numI3eEx = numI3eEx + 1
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
			if (name == 'III4e-extrap'):				
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				exIso.append('III4e-Ex' + str(numIII4eEx))
				exIso.append(xlist2)
				exIso.append(ylist2)
				numIII4eEx = numIII4eEx + 1
			if (name == 'III4e-extrap-full'):				
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				exIsoFull.append('III4e-Ex-full' + str(numIII4eExFull))
				exIsoFull.append(xlist2)
				exIsoFull.append(ylist2)
				numIII4eExFull = numIII4eExFull + 1
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
			if (name == 'V4e-extrap'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				exIso.append('V4e-Ex' + str(numV4eEx))
				exIso.append(xlist2)
				exIso.append(ylist2)
				numV4eEx = numV4eEx + 1
			if (name == 'V4e-extrap-full'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				exIsoFull.append('V4e-Ex-full' + str(numV4eExFull))
				exIsoFull.append(xlist2)
				exIsoFull.append(ylist2)
				numV4eExFull = numV4eExFull + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name =='V4e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				openIso.append('V4e-open' + str(numV4eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				numV4eOpen = numV4eOpen + 1

	return(isopters, scotomas, openIso, exIso, exIsoFull)

def plotField(fileName,scale,isopters,scotomas,openIso,extrapIso,extrapIsoFull,mX,mY,lX,lY,tX,tY,rX,bX,bY,eye,isoptersTested, isoptersInferred, isopter_labels_dict, intersects_dict, root_list, parent_dict):
	print("In plotField: mx=", mX, " my=", mY, " bx=", bX, " by=", bY, " tx=", tX, " ty=", tY, " lx=", lX, " rx=", rX)
	col = color.cmyk(.6,.6,.6,.5)
	text.set(mode="latex")
	#text.set(cls=SingleRunner)
	text.preamble(r"\renewcommand{\familydefault}{\sfdefault}")
	text.preamble(r"\usepackage{color}")
	text.preamble(r"\definecolor{COL}{cmyk}{%g,%g,%g,%g}" % (col.c, col.m, col.y, col.k))
	c = canvas.canvas()

	inferIsopters = 0
	if(len(isoptersTested) == 0 and len(isoptersInferred) > 0):
		inferIsopters = 1
	## DEBUG
	#print(isoptersTested)
	#print(isoptersInferred)
	## if inferred and none annotated as tested, copy the list over
	if(inferIsopters == 1):
		isoptersTested = isoptersInferred.copy()
	else:
		## if any isopters in isoptersInferred that are not in isoptersTested, add those to isoptersTested
		if(len(isoptersInferred) > 0):
			testedDict = {}
			for ii in range(len(isoptersTested)):
				testedDict[isoptersTested[ii]] = 1
			for ii in range(len(isoptersInferred)):
				if isoptersInferred[ii] not in testedDict:
					# isoptersTested.append(isoptersInferred[ii])
					## Throw an error - if they've been annotated they MUST agree
					## likely QA error. Need to be reviewed again.
					print("ERROR: Annotations found for isopters not tested in " + fileName + " " + isoptersInferred[ii])
					sys.exit()

	mX = mX/scale
	mY = mY/scale
	lX = lX/scale
	lY = lY/scale
	bX = bX/scale
	bY = bY/scale
	rX = rX/scale
	deg_length = bY - mY
	d_x = bX - mX
	d_y = bY - mY
	deg_length = math.sqrt((d_x*d_x)+(d_y*d_y))
	chunk = float(deg_length/7)
	## TODD - I think chunk needs to be negative
	if(chunk > 0):
		chunk = chunk * -1.0
	ringSC = []

	## This adjustment is to compenstate the positioning on the Eye label and 
	## isopter legend when the image is upside down
	#print("DEBUG(mX,mY,lX,bY,rX,chunk):", mX, mY, lX, bY, rX, chunk);
	#if(lX > rX):
	#	chunk = chunk * -1

	p = path.line(mX, (mY+chunk*9), mX, (mY-chunk*8.6))
	c.stroke(p, [style.linewidth(.1), color.rgb(1,1,1)])
	p = path.line(mX+chunk*10, mY, mX-chunk*10, mY)
	c.stroke(p, [style.linewidth(.1), color.rgb(1,1,1)])


	working = root_list
	to_trace = ''

	if (len(working) > 0):
		to_trace = working.pop(0)
	else:
		to_trace = ''


	while (to_trace != ''):
		#render to_trace
		x_points = []
		y_points = []
		points = isopter_labels_dict[to_trace]['polygon'].exterior.coords
		for point in points:
			x_points.append(point[0])
			y_points.append(point[1])
		p=path.line((x_points[0]/scale),(y_points[0]/scale),(x_points[1]/scale),(y_points[1]/scale))
		for j in range (2,len(x_points)-1):
			p.append(path.lineto((x_points[j]/scale),(y_points[j]/scale)))
		p.append(path.closepath())
		if (isopter_labels_dict[to_trace]['type'] == 'sc'):
			c.fill(p, [deformer.smoothed(10.0), style.linewidth(.1), isopter_labels_dict[to_trace]['color']])
		if (isopter_labels_dict[to_trace]['type'] == 'iso'):
			sc_count = 0
			parent = parent_dict[to_trace]
			if (parent == {}):
				#print("Filling: ", isopter_labels_dict[to_trace]['name'])
				c.fill(p, [deformer.smoothed(10.0), style.linewidth(.1), color.rgb.white])
				c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), isopter_labels_dict[to_trace]['color']])
			else:
				type_list = []
				flipped_list = []
				type_list.append(isopter_labels_dict[to_trace]['name'])
				while (parent != {}):
					type_list.append(isopter_labels_dict[parent]['name'])
					new_parent = parent_dict[parent]
					parent = new_parent
				flipped_list = type_list[::-1]
				#Make sc list
				sc_color_name = ''
				for i in range(0,len(flipped_list)):
					if ("sc" in flipped_list[i]):
						name = flipped_list[i]
						type = name.split("-")
						level = type[0]
						test = 0
						for j in range(i+1,len(flipped_list)):
							name2 = flipped_list[j]
							type2 = name2.split("-")
							level2 = type2[0]
							if (level == level2):
								test = 1
						if (test == 0):
							sc_color_name = flipped_list[i]
				#print (flipped_list)
				for each in flipped_list:
					if ('sc' in each):
						sc_count = sc_count + 1
					else:
						sc_count = sc_count - 1
						if (sc_count < 0):
							sc_count = 0
				#print (sc_count)
				if (sc_color_name != ''):
					for each in isopter_labels_dict:
						if (isopter_labels_dict[each]['name'] == sc_color_name):
							c.fill(p, [deformer.smoothed(10.0), style.linewidth(.1), isopter_labels_dict[each]['color']])
							c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), isopter_labels_dict[to_trace]['color']])
				elif (sc_count == 0):
					#print("Filling: ", isopter_labels_dict[to_trace]['name'])
					c.fill(p, [deformer.smoothed(10.0), style.linewidth(.1), color.rgb.white])
					c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), isopter_labels_dict[to_trace]['color']])

		for each in intersects_dict[to_trace]["isopters"]:
			working.append(each)
		for each in intersects_dict[to_trace]["scotomas"]:
			working.append(each)
		if (len(working) != 0):
			to_trace = working.pop(0)
		else:
			to_trace = ''

	#Print open isopters
	for i in range(0,len(openIso),3):
		if ('I1e' in openIso[i]):
			p=path.line((openIso[i+1][0]/scale),(openIso[i+2][0]/scale),(openIso[i+1][1]/scale),(openIso[i+2][1]/scale))
			for j in range (2,len(openIso[i+1])):
				p.append(path.lineto((openIso[i+1][j]/scale),(openIso[i+2][j]/scale)))
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), color.rgb.green])
		if ('I2e' in openIso[i]):
			p=path.line((openIso[i+1][0]/scale),(openIso[i+2][0]/scale),(openIso[i+1][1]/scale),(openIso[i+2][1]/scale))
			for j in range (2,len(openIso[i+1])):
				p.append(path.lineto((openIso[i+1][j]/scale),(openIso[i+2][j]/scale)))
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), color.rgb.red])
		if ('I3e' in openIso[i]):
			p=path.line((openIso[i+1][0]/scale),(openIso[i+2][0]/scale),(openIso[i+1][1]/scale),(openIso[i+2][1]/scale))
			for j in range (2,len(openIso[i+1])):
				p.append(path.lineto((openIso[i+1][j]/scale),(openIso[i+2][j]/scale)))
				print (openIso[i+1][j], openIso[i+2][j])
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), color.rgb(.4,.4,.4)])
		if ('I4e' in openIso[i]):
			p=path.line((openIso[i+1][0]/scale),(openIso[i+2][0]/scale),(openIso[i+1][1]/scale),(openIso[i+2][1]/scale))
			for j in range (2,len(openIso[i+1])):
				p.append(path.lineto((openIso[i+1][j]/scale),(openIso[i+2][j]/scale)))
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), color.rgb.blue])
		if ('III4e' in openIso[i]):
			p=path.line((openIso[i+1][0]/scale),(openIso[i+2][0]/scale),(openIso[i+1][1]/scale),(openIso[i+2][1]/scale))
			for j in range (2,len(openIso[i+1])):
				p.append(path.lineto((openIso[i+1][j]/scale),(openIso[i+2][j]/scale)))
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), color.rgb(.8,.4,.1)])
		if ('V4e' in openIso[i]):
			p=path.line((openIso[i+1][0]/scale),(openIso[i+2][0]/scale),(openIso[i+1][1]/scale),(openIso[i+2][1]/scale))
			for j in range (2,len(openIso[i+1])):
				p.append(path.lineto((openIso[i+1][j]/scale),(openIso[i+2][j]/scale)))
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), color.rgb(1,.3,1)])

	#print (extrapIso)
	for i in range(0,len(extrapIso),3):
		if ('I3e' in extrapIso[i]):
			p=path.line((extrapIso[i+1][0]/scale),(extrapIso[i+2][0]/scale),(extrapIso[i+1][1]/scale),(extrapIso[i+2][1]/scale))
			for j in range (2,len(extrapIso[i+1])):
				p.append(path.lineto((extrapIso[i+1][j]/scale),(extrapIso[i+2][j]/scale)))
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), style.linestyle.dashed, color.rgb(.4,.4,.4)])
		if ('III4e' in extrapIso[i]):
			p=path.line((extrapIso[i+1][0]/scale),(extrapIso[i+2][0]/scale),(extrapIso[i+1][1]/scale),(extrapIso[i+2][1]/scale))
			for j in range (2,len(extrapIso[i+1])):
				p.append(path.lineto((extrapIso[i+1][j]/scale),(extrapIso[i+2][j]/scale)))
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), style.linestyle.dashed, color.rgb(.8,.4,.1)])
		if ('V4e' in extrapIso[i]):
			print ("hereeeeeeee")
			p=path.line((extrapIso[i+1][0]/scale),(extrapIso[i+2][0]/scale),(extrapIso[i+1][1]/scale),(extrapIso[i+2][1]/scale))
			for j in range (2,len(extrapIso[i+1])):
				p.append(path.lineto((extrapIso[i+1][j]/scale),(extrapIso[i+2][j]/scale)))
			for k in range(0,len(openIso),3):
				if ('V4e' in openIso[k]):
					print (openIso[k+1])
					p.append(path.lineto((extrapIso[k+1][-1]/scale),(extrapIso[k+2][-1]/scale)))
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), style.linestyle.dashed, color.rgb(1,.3,1)])

	for i in range(0,len(extrapIsoFull),3):
		if ('I3e' in extrapIsoFull[i]):
			p=path.line((extrapIsoFull[i+1][0]/scale),(extrapIsoFull[i+2][0]/scale),(extrapIsoFull[i+1][1]/scale),(extrapIsoFull[i+2][1]/scale))
			for j in range (2,len(extrapIsoFull[i+1])):
				p.append(path.lineto((extrapIsoFull[i+1][j]/scale),(extrapIsoFull[i+2][j]/scale)))
			p.append(path.closepath())
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), style.linestyle.dashed, color.rgb(.4,.4,.4)])
		if ('III4e' in extrapIsoFull[i]):
			p=path.line((extrapIsoFull[i+1][0]/scale),(extrapIsoFull[i+2][0]/scale),(extrapIsoFull[i+1][1]/scale),(extrapIsoFull[i+2][1]/scale))
			for j in range (2,len(extrapIsoFull[i+1])):
				p.append(path.lineto((extrapIsoFull[i+1][j]/scale),(extrapIsoFull[i+2][j]/scale)))
			p.append(path.closepath())
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), style.linestyle.dashed, color.rgb(.8,.4,.1)])
		if ('V4e' in extrapIsoFull[i]):
			p=path.line((extrapIsoFull[i+1][0]/scale),(extrapIsoFull[i+2][0]/scale),(extrapIsoFull[i+1][1]/scale),(extrapIsoFull[i+2][1]/scale))
			for j in range (2,len(extrapIsoFull[i+1])):
				p.append(path.lineto((extrapIsoFull[i+1][j]/scale),(extrapIsoFull[i+2][j]/scale)))
			p.append(path.closepath())
			c.stroke(p, [deformer.smoothed(10.0), style.linewidth(.2), style.linestyle.dashed, color.rgb(1,.3,1)])
		
	#Draw circles
	# mX,mY,lX,tY,rX
	#radius = lX - mX
	r_x = lX - mX
	r_y = lY - mY
	radius = math.sqrt((r_x*r_x) + (r_y*r_y))
	#print("DEBUG: radius = ",radius,"  chunk = ",chunk)
	#Tick marks
	deg_length = bY - mY
	for i in range(1,8):
		radius = i*chunk
		if i == 3:
			c.stroke(path.circle(mX, mY, radius),[ style.linewidth(.1), color.transparency(.45)])
		else:
			c.stroke(path.circle(mX, mY, radius),[ style.linewidth(.1), color.transparency(.75)])
	c.stroke(path.circle(mX, mY, chunk*.5),[ style.linewidth(.1), color.transparency(.75)])

	#Bottom distances
	distance1 = math.tan(math.radians(15))*(7.4*chunk)
	distance2 = math.tan(math.radians(30))*(7.4*chunk)
	#Top distances
	distance3 = math.tan(math.radians(165))*(7*chunk)
	distance4 = math.tan(math.radians(150))*(7*chunk)
	#Radial Lines
	x10_15 = (.5*chunk*(math.sin(math.radians(15)))) + mX
	y10_15 = (.5*chunk*(math.cos(math.radians(15)))) + mY
	y90_15 = (7.4*chunk) + mY
	x90_15 = mX + distance1
	p = path.line(x10_15, y10_15, x90_15, y90_15)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_30 = (.5*chunk*(math.sin(math.radians(30)))) + mX
	y10_30 = (.5*chunk*(math.cos(math.radians(30)))) + mY
	x90_30 = mX + distance2
	y90_30 = (7.4*chunk) + mY
	p = path.line(x10_30, y10_30, x90_30, y90_15)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_45 = (.5*chunk*(math.sin(math.radians(45)))) + mX
	y10_45 = (.5*chunk*(math.cos(math.radians(45)))) + mY
	x90_45 = (9*chunk*(math.sin(math.radians(45)))) + mX
	y90_45 = (9*chunk*(math.cos(math.radians(45)))) + mY
	p = path.line(x10_45, y10_45, x90_45, y90_45)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_60 = (.5*chunk*(math.sin(math.radians(60)))) + mX
	y10_60 = (.5*chunk*(math.cos(math.radians(60)))) + mY
	x90_60 = (9*chunk*(math.sin(math.radians(60)))) + mX
	y90_60 = (9*chunk*(math.cos(math.radians(60)))) + mY
	p = path.line(x10_60, y10_60, x90_60, y90_60)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_75 = (.5*chunk*(math.sin(math.radians(75)))) + mX
	y10_75 = (.5*chunk*(math.cos(math.radians(75)))) + mY
	x90_75 = (9*chunk*(math.sin(math.radians(75)))) + mX
	y90_75 = (9*chunk*(math.cos(math.radians(75)))) + mY
	p = path.line(x10_75, y10_75, x90_75, y90_75)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_105 = (.5*chunk*(math.sin(math.radians(105)))) + mX
	y10_105 = (.5*chunk*(math.cos(math.radians(105)))) + mY
	x90_105 = (9*chunk*(math.sin(math.radians(105)))) + mX
	y90_105 = (9*chunk*(math.cos(math.radians(105)))) + mY
	p = path.line(x10_105, y10_105, x90_105, y90_105)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_120 = (.5*chunk*(math.sin(math.radians(120)))) + mX
	y10_120 = (.5*chunk*(math.cos(math.radians(120)))) + mY
	x90_120 = (9*chunk*(math.sin(math.radians(120)))) + mX
	y90_120 = (9*chunk*(math.cos(math.radians(120)))) + mY
	p = path.line(x10_120, y10_120, x90_120, y90_120)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_135 = (.5*chunk*(math.sin(math.radians(135)))) + mX
	y10_135 = (.5*chunk*(math.cos(math.radians(135)))) + mY
	x90_135 = (9*chunk*(math.sin(math.radians(135)))) + mX
	y90_135 = (9*chunk*(math.cos(math.radians(135)))) + mY
	p = path.line(x10_135, y10_135, x90_135, y90_135)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_150 = (.5*chunk*(math.sin(math.radians(150)))) + mX
	y10_150 = (.5*chunk*(math.cos(math.radians(150)))) + mY
	y90_150 = mY - (7*chunk)
	x90_150 = mX - distance4
	p = path.line(x10_150, y10_150, x90_150, y90_150)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_165 = (.5*chunk*(math.sin(math.radians(165)))) + mX
	y10_165 = (.5*chunk*(math.cos(math.radians(165)))) + mY
	#x90_165 = (9*chunk*(math.sin(math.radians(165)))) + mX
	#y90_165 = (9*chunk*(math.cos(math.radians(165)))) + mY
	y90_165 = mY - (7*chunk)
	x90_165 = mX - distance3
	p = path.line(x10_165, y10_165, x90_165, y90_165)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_195 = (.5*chunk*(math.sin(math.radians(195)))) + mX
	y10_195 = (.5*chunk*(math.cos(math.radians(195)))) + mY
	y90_195 = mY - (7*chunk)
	x90_195 = mX + distance3
	p = path.line(x10_195, y10_195, x90_195, y90_195)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_210 = (.5*chunk*(math.sin(math.radians(210)))) + mX
	y10_210 = (.5*chunk*(math.cos(math.radians(210)))) + mY
	#x90_210 = (9*chunk*(math.sin(math.radians(210)))) + mX
	#y90_210 = (9*chunk*(math.cos(math.radians(210)))) + mY
	y90_210 = mY - (7*chunk)
	x90_210 = mX + distance4
	p = path.line(x10_210, y10_210, x90_210, y90_210)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_225 = (.5*chunk*(math.sin(math.radians(225)))) + mX
	y10_225 = (.5*chunk*(math.cos(math.radians(225)))) + mY
	x90_225 = (9*chunk*(math.sin(math.radians(225)))) + mX
	y90_225 = (9*chunk*(math.cos(math.radians(225)))) + mY
	p = path.line(x10_225, y10_225, x90_225, y90_225)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_240 = (.5*chunk*(math.sin(math.radians(240)))) + mX
	y10_240 = (.5*chunk*(math.cos(math.radians(240)))) + mY
	x90_240 = (9*chunk*(math.sin(math.radians(240)))) + mX
	y90_240 = (9*chunk*(math.cos(math.radians(240)))) + mY
	p = path.line(x10_240, y10_240, x90_240, y90_240)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_255 = (.5*chunk*(math.sin(math.radians(255)))) + mX
	y10_255 = (.5*chunk*(math.cos(math.radians(255)))) + mY
	x90_255 = (9*chunk*(math.sin(math.radians(255)))) + mX
	y90_255 = (9*chunk*(math.cos(math.radians(255)))) + mY
	p = path.line(x10_255, y10_255, x90_255, y90_255)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_285 = (.5*chunk*(math.sin(math.radians(285)))) + mX
	y10_285 = (.5*chunk*(math.cos(math.radians(285)))) + mY
	x90_285 = (9*chunk*(math.sin(math.radians(285)))) + mX
	y90_285 = (9*chunk*(math.cos(math.radians(285)))) + mY
	p = path.line(x10_285, y10_285, x90_285, y90_285)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_300 = (.5*chunk*(math.sin(math.radians(300)))) + mX
	y10_300 = (.5*chunk*(math.cos(math.radians(300)))) + mY
	x90_300 = (9*chunk*(math.sin(math.radians(300)))) + mX
	y90_300 = (9*chunk*(math.cos(math.radians(300)))) + mY
	p = path.line(x10_300, y10_300, x90_300, y90_300)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_315 = (.5*chunk*(math.sin(math.radians(315)))) + mX
	y10_315 = (.5*chunk*(math.cos(math.radians(315)))) + mY
	x90_315 = (9*chunk*(math.sin(math.radians(315)))) + mX
	y90_315 = (9*chunk*(math.cos(math.radians(315)))) + mY
	p = path.line(x10_315, y10_315, x90_315, y90_315)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_330 = (.5*chunk*(math.sin(math.radians(330)))) + mX
	y10_330 = (.5*chunk*(math.cos(math.radians(330)))) + mY
	x90_330 = mX - distance2
	y90_330 = (7.4*chunk) + mY
	p = path.line(x10_330, y10_330, x90_330, y90_330)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x10_345 = (.5*chunk*(math.sin(math.radians(345)))) + mX
	y10_345 = (.5*chunk*(math.cos(math.radians(345)))) + mY
	x90_345 = mX - distance1
	y90_345 = (7.4*chunk) + mY
	p = path.line(x10_345, y10_345, x90_345, y90_345)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	#Cross lines
	#Top and Bottom
	p = path.line(mX, mY, mX, (mY-(chunk*7)))
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	p = path.line(mX, mY, mX, (mY+(chunk*7.4)))
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	#Left and Right
	p = path.line(mX, mY, (mX-(chunk*9)), mY)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	p = path.line(mX, mY, (mX+(chunk*9)), mY)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	#Draw top and bottom lines
	x1 = (9*chunk*(math.sin(math.radians(320.6)))) + mX
	y1 = mY - (9*chunk*(math.cos(math.radians(320.8))))
	x2 = (9*chunk*(math.sin(math.radians(39.2)))) + mX
	y2 = mY - (9*chunk*(math.cos(math.radians(39.2))))
	p = path.line(x1, y1, x2, y2)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	x1 = (9*chunk*(math.sin(math.radians(145.6)))) + mX
	y = mY+(chunk*7.43)
	x2 = (9*chunk*(math.sin(math.radians(214.4)))) + mX
	p = path.line(x1, y, x2, y)
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	#Bottom left arcs
	p = path.path(path.arc(mX,mY,(chunk*9),0,55.6))
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	p = path.path(path.arc(mX,mY,(chunk*8),0,68))
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	#Bottom right arcs
	p = path.path(path.arc(mX,mY,(chunk*9),124.4,180))
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	p = path.path(path.arc(mX,mY,(chunk*8),112,180))
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	#Top right arcs
	p = path.path(path.arc(mX,mY,(chunk*9),180,230.7))
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	p = path.path(path.arc(mX,mY,(chunk*8),180,241))
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	#Top left arcs
	p = path.path(path.arc(mX,mY,(chunk*9),309,360))
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])
	p = path.path(path.arc(mX,mY,(chunk*8),299,360))
	c.stroke(p, [style.linewidth(.1), color.transparency(.75)])

	#Number labels
	c.text(mX, ( mY+chunk-.4), r"\textcolor{COL}{1 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text(mX, ( mY+(chunk*2)-.4), r"\textcolor{COL}{2 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text(mX, ( mY+(chunk*3)-.4), r"\textcolor{COL}{3 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text(mX, ( mY+(chunk*4)-.4), r"\textcolor{COL}{4 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text(mX, ( mY+(chunk*5)-.4), r"\textcolor{COL}{5 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text(mX, ( mY+(chunk*6)-.4), r"\textcolor{COL}{6 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text(mX, ( mY+(chunk*7)-.15), r"\textcolor{COL}{7 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.top])
	c.text(mX, ( mY-chunk+.4), r"\textcolor{COL}{1 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.top])
	c.text(mX, ( mY-(chunk*2)+.4), r"\textcolor{COL}{2 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.top])
	c.text(mX, ( mY-(chunk*3)+.4), r"\textcolor{COL}{3 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.top])
	c.text(mX, ( mY-(chunk*4)+.4), r"\textcolor{COL}{4 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.top])
	c.text(mX, ( mY-(chunk*5)+.4), r"\textcolor{COL}{5 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.top])
	c.text(mX, ( mY-(chunk*6)+.4), r"\textcolor{COL}{6 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.top])
	c.text(mX, ( mY-(chunk*7)-.1), r"\textcolor{COL}{7 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.top])
	c.text((mX+chunk), (mY+.1), r"\textcolor{COL}{1 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX+(chunk*2)), (mY+.1), r"\textcolor{COL}{2 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX+(chunk*3)), (mY+.1), r"\textcolor{COL}{3 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX+(chunk*4)), (mY+.1), r"\textcolor{COL}{4 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX+(chunk*5)), (mY+.1), r"\textcolor{COL}{5 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX+(chunk*6)), (mY+.1), r"\textcolor{COL}{6 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX+(chunk*7)), (mY+.1), r"\textcolor{COL}{7 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX+(chunk*8)), (mY+.1), r"\textcolor{COL}{8 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX+(chunk*9)+.45), (mY+.1), r"\textcolor{COL}{90}", [text.size(1), text.halign.boxright, text.halign.flushright, text.valign.bottom])
	c.text((mX-chunk), (mY+.1), r"\textcolor{COL}{1 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX-(chunk*2)), (mY+.1), r"\textcolor{COL}{2 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX-(chunk*3)), (mY+.1), r"\textcolor{COL}{3 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX-(chunk*4)), (mY+.1), r"\textcolor{COL}{4 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX-(chunk*5)), (mY+.1), r"\textcolor{COL}{5 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX-(chunk*6)), (mY+.1), r"\textcolor{COL}{6 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX-(chunk*7)), (mY+.1), r"\textcolor{COL}{7 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX-(chunk*8)), (mY+.1), r"\textcolor{COL}{8 0}", [text.size(1), text.halign.boxcenter, text.halign.flushcenter, text.valign.bottom])
	c.text((mX-(chunk*9)-.45), (mY+.1), r"\textcolor{COL}{90}", [text.size(1), text.halign.boxleft, text.halign.flushleft, text.valign.bottom])


	#Add text for left vs. right eye
	#x = mX - (9.3*chunk*(math.sin(math.radians(137))))
	#y = mY - (9.5*chunk*(math.cos(math.radians(137))))
	x = mX - (9.3*chunk*(math.sin(math.radians(137)))) - 0.4
	y = mY - (9.5*chunk*(math.cos(math.radians(125)))) - 3.4
	c.text(x,y, r""+eye, [text.size(5), text.halign.boxleft, text.halign.flushleft]) # text.fontmap = PDFstdfont("Helvetica")]) 

        ## print isopters
	## First, determine how many there are so we can allocate the proper amount of height for the legend

	## need to sort the isoptersTested list/array
	isoptersTested.sort()

	num_isopters = len (isoptersTested)
	x = mX - (9.3*chunk*(math.sin(math.radians(125)))) - .4
	boxX1 = x
	## will need a dictionary for colors based upon isopter
	color_dict = { "I1e": color.rgb.green, "I2e": color.rgb.red, "I3e": color.rgb(.4,.4,.4),
		"I4e": color.rgb.blue, "III4e": color.rgb(.8,.5,.2), "V4e": color.rgb(1,.3,1) }
	for ii in range(num_isopters):
		y = mY - (9.5*chunk*(math.cos(math.radians(125)))) - 3.4 + (0.6 * (num_isopters - ii - 1))
		boxY1 = y -.1
		mycolor = color_dict[isoptersTested[ii]]
		c.stroke(path.rect(x,y,.5,.5), [mycolor, deco.filled([mycolor])])
		c.text(x+1,y+.1, r""+isoptersTested[ii], [text.size(2), text.halign.boxleft, text.halign.flushleft, text.valign.bottom]) 

	## draw the bounding box around the legend
	legend_height = 0.1 + (0.6 * num_isopters)
	c.stroke(path.rect(boxX1-.1, boxY1, 2.3, legend_height))

	y = mY - (9.5*chunk*(math.cos(math.radians(125)))) - 4.1
	x = x - 0.15
	#c.text(x,y, r"Isopters tested", [text.size(2), text.halign.boxleft, text.halign.flushleft, text.valign.bottom]) 
	#y = y - 0.6
	if(inferIsopters == 1):
		c.text(x,y, r"*inferred from", [text.size(2), text.halign.boxleft, text.halign.flushleft, text.valign.bottom]) 
		c.text(x+0.2,y-0.6, r"annotated isopters", [text.size(2), text.halign.boxleft, text.halign.flushleft, text.valign.bottom]) 


	c.writePDFfile(fileName)

def construct_hos(isopters, scotomas, openIso, mX, mY, tX, tY, bX, bY, rX, rY, lX, lY):

	alphabet_string = string.ascii_lowercase
	polygons_iso = []
	polygons_sc = []

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
	all_isopters = polygons_iso + polygons_sc

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
			isopter_labels_dict[letter_name]['color'] = color.rgb.green
		if "I2e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
			isopter_labels_dict[letter_name]['color'] = color.rgb.red
		if "I3e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
			isopter_labels_dict[letter_name]['color'] = color.rgb(.4,.4,.4)
		if "I4e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
			isopter_labels_dict[letter_name]['color'] = color.rgb.blue
		if "III4e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
			isopter_labels_dict[letter_name]['color'] = color.rgb(.8,.4,.1)
		if "V4e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
			isopter_labels_dict[letter_name]['color'] =  color.rgb(1,.3,1)


	#Find all intersections
	intersects_dict = {}
	for each in isopter_labels_dict:
		intersect_list_sc = []
		intersect_list_iso = []
		intersects_dict[each] = {}
		intersects_dict[each]["isopters"] = {}
		intersects_dict[each]["scotomas"] = {}
		for other in isopter_labels_dict:
			if (each != other):
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

	#Need to make object that contains the parent of each entry
	parent_dict = {}
	for each in isopter_labels_dict:
		parent_dict[each] = {}
		for other in intersects_dict:
			if each in intersects_dict[other]["isopters"]:
				parent_dict[each] = other
			if each in intersects_dict[other]["scotomas"]:
				parent_dict[each] = other

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

	return (isopter_labels_dict, intersects_dict, root_list, parent_dict)


def parseAnnotation(inputfile, output_dir, output_filepath, use_log):
	inputfile = inputfile
	dom = xml.dom.minidom.parse(inputfile)
	images = dom.getElementsByTagName("annotation-per-image")
	log_fh = None
	logName = "log.txt"

	#if (use_log):
	#	log_fh = open(logName,"w")
	#	log_fh.write("Image" + "\t" + "Sufficient Fiducials" + "\t" + "Marked Curves" + "\t" + "Eye Annotated" + "\t" + "PDF Generated" + "\n")
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
			#print("DEBUG: AveRRatio = ", avgRRatio)
			#Rotate fiducial points
			rotatedFids, pts = rotateFiducials(avgRotation,avgRRatio,mX,mY,tX,tY,rX,rY,bX,bY,lX,lY)
			#Will rotate points and flip y-axis
			isopters, scotomas, openIso, extrapIso, extrapIsoFull = collectIsopters(each,avgRotation,avgRRatio,mX,mY)
			isopter_labels_dict, intersects_dict, root_list, parent_dict = construct_hos(isopters, scotomas, openIso, mX, mY, tX, tY, bX, bY, rX, rY, lX, lY)
			isopters_tested, isopters_inferred, isopters_tested_annotated = getIsoptersTested(each)
			if(len(isopters_tested) == 0 and len(isopters_inferred) == 0):
				## nothing to graph, so get out of here.
				print("ERROR: No isopters tested, and no isopters annotated in XML.")
				break
			if(isopters_tested_annotated == 1 and len(isopters_tested) == 0 and len(isopters_inferred) > 0):
				## throw an error and quit
				print("ERROR: isotpers_tested tag is present, all are marked untested, but isopters found in XML.")
				break
			#plotField(fileName,scale,isopters,scotomas,openIso,mX,mY,lX,lY,tX,tY,rX,eye,isopters_tested, isopters_inferred)
			plotField(fileName,scale,isopters,scotomas,openIso,extrapIso,extrapIsoFull,mX,mY,lX,lY,tX,tY,rX,bX,bY,eye,isopters_tested, isopters_inferred, isopter_labels_dict, intersects_dict, root_list, parent_dict)
			if ((len(isopters) > 1) or (len(scotomas) > 1) or (len(openIso) > 1)):
				markedCurves = "Yes"
			else:
				markedCurves = "No"
			if (eye == "OS" or eye == "OD"):
				eyeAnnotate = "Yes"
			else:
				eyeAnnotate = "No"
			PDFGenerated = "Yes"
			#if (use_log):
			#	log_fh.write(fileName + "\t" + sufficientFid + "\t" + markedCurves + "\t" + eyeAnnotate + "\t" + PDFGenerated + "\n")
		else:
			sufficientFid = "No"
			markedCurves = "NA"
			eyeAnnotate = "NA"
			PDFGenerated = "No"
			#if (use_log):
			#	log_fh.write(fileName + "\t" + sufficientFid + "\t" + markedCurves + "\t" + eyeAnnotate + "\t" + PDFGenerated + "\n")

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
		elif opt == "--output-file-path":
			output_filepath = arg
			#print ('Output filepath is ', output_filepath)
		elif opt in ("-o", "--output-dir"):
			output_dir = arg
			#print ('Output Directory is ', output_dir)
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
