#!/usr/bin/python

import math
import os, sys, getopt, pprint, re, math
import xml.dom.minidom
import numpy as np
import os
import datetime
import dateutil
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import string
import matplotlib.pyplot as plt

#from pprint import pprint
#from pyx import *
from shapely.geometry import Polygon
from shapely.geometry import Point
from shapely.geometry import LinearRing
from shapely.geometry import LineString
from shapely.geometry import GeometryCollection
from shapely.geometry import MultiPoint
from shapely.geometry import MultiLineString
from shapely import geometry
from skimage.morphology import skeletonize
from skimage import data
from skimage.util import invert


## find three closest children
def find_three_closest_children(isopter_labels_dict, below_iso, all_children, test_point):
	#print("DEBUG: Finding closest children of ", isopter_labels_dict[below_iso]['name'])

	child_dict = {}
	## find distance from test_point to below_iso
	tt_point = (test_point.x, test_point.y)
	distance_below, point_below = find_distance_to_isopter(isopter_labels_dict, tt_point, below_iso)
	## find distance to each child, retain a list of distances for all children
	min_dist_to_child = 100000
	min_dist_to_child_iso = ''
	for each in all_children:
		distance_each, point_each = find_distance_to_isopter(isopter_labels_dict, tt_point, each)
		child_dict[each] = {}
		child_dict[each]['dist'] = distance_each
		child_dict[each]['point'] = point_each



## find two closest children
## need a function that determines if there are two children isopters that are closest to the test point,
## and for which a line from the test point to the child does not pass through any other children
## that function should return the IDs of the two children if two are found.
## Otherwise, return the ID of the closest child as child 1, and empty string as child2
def find_two_closest_children(isopter_labels_dict, below_iso, all_children, test_point):
	#print("DEBUG: Finding closest children of ", isopter_labels_dict[below_iso]['name'])

	child_dict = {}
	## find distance from test_point to below_iso
	tt_point = (test_point.x, test_point.y)
	distance_below, point_below = find_distance_to_isopter(isopter_labels_dict, tt_point, below_iso)
	## find distance to each child, retain a list of distances for all children
	min_dist_to_child = 100000
	min_dist_to_child_iso = ''
	for each in all_children:
		distance_each, point_each = find_distance_to_isopter(isopter_labels_dict, tt_point, each)
		child_dict[each] = {}
		child_dict[each]['dist'] = distance_each
		child_dict[each]['point'] = point_each

	## now, find the two children that are the closest (smallest distance) and which have an un-occluded "view" of the point

	child_list_by_distance = sorted(child_dict, key=lambda x: child_dict[x]['dist'])
	c1 = child_list_by_distance.pop(0)
	c2 = ''
	## First, test min_1 vs min_2
	x1 = test_point.x
	y1 = test_point.y
	for i in child_list_by_distance:
		## check if line from pt to i goes through c1
		p2 = child_dict[i]['point']
		x2 = p2[0]
		y2 = p2[1]
		line = LineString([(x1,y1),(x2,y2)])
		c1_polygon = isopter_labels_dict[c1]['polygon']
		## TODD - test with values other than 0 (0.5, 1.0, 2.0)
		c1_polygon = c1_polygon.buffer(0)
		inter = c1_polygon.intersection(line)
		#print("DEBUG: INTER =", inter)
		## FLAG = 1 --> INTERSECTION; FLAG = 0 --> NO INTERSECTION
		flag = 1	
		if (isinstance(inter,LineString)):
			#print("DEBUG: LineString len(INTER.COORDS) =", len(inter.coords))
			if(len(inter.coords) == 0):
				flag = 0
		elif (isinstance(inter,GeometryCollection)):
			#print("DEBUG: GeomColl len(INTER.GEOM) =", len(inter.geoms))
			flag = 0
		elif (isinstance(inter,Point)):
			#print("DEBUG: Point len(INTER.COORDS) =", len(inter.coords))
			flag = 0
		if(flag == 0):
			#print("DEBUG: NO INTERSECTION")
			#print("C1 = %s, I = %s" % (c1,i))
			c2 = i
			return(c1, c2)
	## Fall-through. Only the closest child has an unoccluded view of the test point
	return(c1, '')




def construct_hos(isopters, scotomas, circle_points, circle_degs, mX, mY, tX, tY, bX, bY, rX, rY, lX, lY, degree_len, output_file, z_scale, eye):

	alphabet_string = string.ascii_lowercase
	polygons_iso = []
	polygons_sc = []
	
	#isopters and scotomas to Polygons
	for i in range(0,len(isopters),4):
		point_list = []
		for j in range(0,len(isopters[i+1])):
			point_list.append((isopters[i+1][j], isopters[i+2][j]))
		polygon = Polygon(point_list)
		polygons_iso.append(isopters[i])
		polygons_iso.append(polygon)
		polygons_iso.append(isopters[i+3][0])
	for i in range(0,len(scotomas),4):
		point_list = []
		for j in range(0,len(scotomas[i+1])):
			point_list.append((scotomas[i+1][j], scotomas[i+2][j]))
		scotoma = Polygon(point_list)
		polygons_sc.append(scotomas[i])
		polygons_sc.append(scotoma)
		polygons_sc.append(scotomas[i+3][0])
	all_isopters = polygons_iso + polygons_sc

	isopter_labels_dict = {}
	letters_count = 0
	for i in range(0,len(all_isopters),3):
		letter_name = alphabet_string[letters_count]
		letters_count = letters_count + 1
		isopter_labels_dict[letter_name] = {}
		isopter_labels_dict[letter_name]['name'] = all_isopters[i]
		isopter_labels_dict[letter_name]['polygon'] = all_isopters[i+1]
		isopter_labels_dict[letter_name]['height'] = all_isopters[i+2]
		if "I1e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
		if "I2e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
		if "I3e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
		if "I4e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
		if "III4e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
		if "V4e" in isopter_labels_dict[letter_name]['name']:
			isopter_labels_dict[letter_name]['type'] = "iso"
			if "sc" in isopter_labels_dict[letter_name]['name']:
				isopter_labels_dict[letter_name]['type'] = "sc"
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

	#print("DEBUG: Before pruning: ID = ", intersects_dict)
	#Prune HOS
#	for each in intersects_dict:
#		sc_to_del = []
#		for each_isopter in intersects_dict[each]["isopters"]:
#			if len(intersects_dict[each_isopter]["scotomas"]) > 0:
#				sc_to_del.append(intersects_dict[each_isopter]["scotomas"][0])
#		if (len(sc_to_del)) >= 1:
#			for each_sc in sc_to_del:
#				if (each_sc in intersects_dict[each]["scotomas"]):
#					intersects_dict[each]["scotomas"].remove(each_sc)
#	for each in intersects_dict:
#		iso_to_del = []
#		for each_isopter in intersects_dict[each]["isopters"]:
#			for other_isopter in intersects_dict[each_isopter]["isopters"]:
#				if (other_isopter in intersects_dict[each]["isopters"]):
#					iso_to_del.append(other_isopter)
#		for each_iso in iso_to_del:
#			if (each_iso in intersects_dict[each]["isopters"]):
#				intersects_dict[each]["isopters"].remove(each_iso)

	for each in intersects_dict:
		## get a list of all children of each
		icmain, scmain = find_children(isopter_labels_dict, intersects_dict, each)
		allmain = icmain + scmain
		for isopter_child in allmain:
			## get all children of isopter_child
			ic, sc = find_children(isopter_labels_dict, intersects_dict, isopter_child)
			all = ic + sc
			## if the child exist in isopter_child and each, remove them them from each
			for test in all:
				if test in intersects_dict[each]["isopters"]:
					intersects_dict[each]["isopters"].remove(test)
				if test in intersects_dict[each]["scotomas"]:
					intersects_dict[each]["scotomas"].remove(test)



	
	#print("DEBUG: After pruning: ID = ", intersects_dict)

	#Need to find isopters to make skeletons of 
	#print("DEBUG: Finding isopters for which skeletons are needed.")
	for each in isopter_labels_dict:
		#print("DEBUG: each = %s" % each)
		isopter_children, scotoma_children = find_children(isopter_labels_dict, intersects_dict, each)
		## Let's try making skeletons for everyone! (TODD)
		skel_points_x, skel_points_y, a, skeleton = make_skeletons(isopter_labels_dict, each, scotoma_children, isopter_children,  lX, rX, tY, bY)
		#display_skeleton(a,skeleton)
		skel_height = skeleton_z(isopter_labels_dict, z_scale,  each)
		skel_points_z = [skel_height] * len(skel_points_x)
		isopter_labels_dict[each]['skeleton'] = (skel_points_x, skel_points_y, skel_points_z)
		## END TODD- PREVIOUSLY
		#if len(isopter_children) == 0 and len(scotoma_children) == 0:
		#	skel_points_x, skel_points_y, a, skeleton = make_skeletons(isopter_labels_dict, each, scotoma_children, isopter_children, lX, rX, tY, bY)
		#	skel_height = skeleton_z(isopter_labels_dict, z_scale, each)
		#	#print (skel_height)
		#	skel_points_z = [skel_height] * len(skel_points_x)
		#	isopter_labels_dict[each]['skeleton'] = (skel_points_x, skel_points_y, skel_points_z)
		#if len(scotoma_children) > 0:
		#	skel_points_x, skel_points_y, a, skeleton = make_skeletons(isopter_labels_dict, each, scotoma_children, isopter_children,  lX, rX, tY, bY)
		#	skel_height = skeleton_z(isopter_labels_dict, z_scale,  each)
		#	skel_points_z = [skel_height] * len(skel_points_x)
		#	isopter_labels_dict[each]['skeleton'] = (skel_points_x, skel_points_y, skel_points_z)
		#	#print("HERE 3:")
		#	#print(isopter_labels_dict[each])
		#	#print(isopter_labels_dict[each]['polygon'].exterior.coords.xy)
	#print("DEBUG: Done making skeletons")

	#print ("isopter_labels_dict")
	#print (isopter_labels_dict)
#	for each in isopter_labels_dict:
#		print (isopter_labels_dict[each])
#		#print (isopter_labels_dict[each]['type'])
#		points = isopter_labels_dict[each]['polygon'].exterior.coords
#		for point in points:
#			print (point[0], point[1])
#		#print (isopter_labels_dict[each]['polygon'])

	#print (isopter_labels_dict)
	#open_len = len(openIso)

	#print (circle_points)
	circle_points_list = []
	for each_circle in circle_points:
		for each_point in each_circle:
			#print(each_point)
			circle_points_list.append(each_point)
	#print (circle_degs)
	circle_degs_list = []
	for each_deg in circle_degs:
		for each_deg2 in each_deg:
			#print(each_point)
			circle_degs_list.append(each_deg2)

#	Testing if all points are here
#	for each in circle_points_list:
#		print (each)

	#print ("\n")
	#print (circle_points[9])
	#print (circle_points[10])
	#circle_points_height = [0] * len(circle_points_list)


	#Make root list
	root_list = []
	for each in isopter_labels_dict:
		root_list.append(each)

	#print (root_list)

	for each in isopter_labels_dict:
		for other in intersects_dict:
			if (each in intersects_dict[other]['isopters']):
				if (each in root_list):
					root_list.remove(each)
			if (each in intersects_dict[other]['scotomas']):
				if (each in root_list):
					root_list.remove(each)


	################################################
	#Need to test to calculate area before moving on
	total_area = 0
	for i in range(0,90,1):
		for j in range(0,360,1):
			n = math.radians(i)
			s = math.radians(i+1)
			e = math.radians(j)
			w = math.radians(j+1)
			area = (math.cos(s)-math.cos(n))*(e-w)
			total_area = total_area + area
#	print (total_area)

	###############################################			
	#Interpolate for each point
	areas = []
	heights = []
	volumes = []
	points = []  ## store X-Y coordinates for each point
	isopter_area = { "V4e" : 0, "III4e": 0, "I4e": 0, "I3e": 0, "I2e": 0, "I1e": 0}
	isopter_height = {"V4e": 0.5, "III4e": 1.5, "I4e": 2.5, "I3e": 3.0, "I2e": 3.5, "I1e": 4.0 }
	total_volume = 0
#	for i in range(40,42,2):
#		for j in range(180,181,1):
	for i in range(0,100,2):
		for j in range(0,360,2):
			#print("DEBUG: r = %d, theta = %0d" % (i, j))
			n = math.radians(i)
			s = math.radians(i+2)
			e = math.radians(j)
			w = math.radians(j+2)
			area = (math.cos(s)-math.cos(n))*(e-w)
			above_iso = ''
			below_iso = ''
			point_height = 0
			## Original, but needs flip on Y axis
			circle_point = ((math.cos(e)*(i*degree_len))+mX),((math.sin(e)*(i*degree_len))+mY)
			## Nope. Get the same result. Need to flip the Y later, pre-graphing
			#circle_point = ((math.cos(e)*(i*degree_len))+mX),mY-((math.sin(e)*(i*degree_len)))
			test_point = Point(circle_point[0], circle_point[1])
			#print("DEBUG: x = %0.3f, y = %0.3f" % (circle_point[0], circle_point[1]))
			for isopter in root_list:
				polygon = isopter_labels_dict[isopter]['polygon']
				#print("DEBUG: polygon = ", polygon)
				point_test = test_point.within(polygon)
				if (point_test == True):
					#print("DEBUG: inside_root: ", isopter_labels_dict[isopter]['name'])
					bound_isopter = isopter_labels_dict[isopter]	
					bounding_isopters = find_boundary(isopter_labels_dict, intersects_dict, bound_isopter, circle_point)
					#print("DEBUG: after find_boundary")
					#Going to need to do something if the returned isopter is a scotoma?
					below_iso = bounding_isopters[0]
					below_iso_name = below_iso['name']
					for letter in isopter_labels_dict:
						letter_name = isopter_labels_dict[letter]['name']
						if (below_iso_name == letter_name):
							below_iso = letter
					above_iso = bounding_isopters[1]

			## above_iso should be the isopter child of below_iso that is closest to test point. But we need
			## more information, if cases where there are multiple children. (It's complicated)

			## if below_iso = '' then not in an isopter, set the height to 0
			if below_iso == '':
				#print("DEBUG: Not in a root isopter. Z = 0")
				point_height = 0

			##Going to need to extrapolate up if you have an above isopter
			else:
				## compute the children of the below_iso
				#print("DEBUG: Before finding children of ", isopter_labels_dict[below_iso]['name'])
				isopter_children, scotoma_children = find_children(isopter_labels_dict, intersects_dict, below_iso)
				all_children = isopter_children + scotoma_children

				#print("PRE_CASE: below_iso = %s, above_iso = %s" % (below_iso, above_iso))
				## if below_iso has more than one child...
				if(len(all_children) > 1):
				#if(len(all_children) > 1 or len(scotoma_children) >0):
					#print("CASE1")
					## convert skeleton into a MultiPoint and determine which child (or below_iso) has an unoccluded view to the test point
					## if none, set height as the skeleton height
					## else, interpolate_two_skel()

					x_points = isopter_labels_dict[below_iso]['skeleton'][0]
					y_points = isopter_labels_dict[below_iso]['skeleton'][1]
					sk_points = []
					for ii in range(0,len(x_points)):
						p = Point(x_points[ii], y_points[ii])
						sk_points.append(p)
					## create the MultiPoint object
					skeleton_mp = MultiPoint(sk_points)
					#print("DEBUG: MP_SKELETON=", skeleton_mp)
					skeleton_mp = skeleton_mp.buffer(1.0)
					
					iso_for_height_check = ''
					distance_for_height_check = 100000

					## add below_iso onto all_children list
					all_children = all_children + [below_iso]

					distance_skel, point_skel = find_distance_to_skeleton(isopter_labels_dict, circle_point, below_iso)
					for iso in all_children:
						## check below_iso, using find_distance to get the closest point
						distance_iso, point_iso = find_distance_to_isopter(isopter_labels_dict, circle_point, iso)
						## create linestring from circle_point to point_below
						iso_ls = LineString([(circle_point[0], circle_point[1]), (point_iso[0], point_iso[1])])
						inter = iso_ls.intersection(skeleton_mp)

						flag = 0   ## no intersection yet
						if(isinstance(inter,LineString)):
							if(len(inter.coords) > 0):
								flag = 1
						elif(isinstance(inter,GeometryCollection)):
							if(len(inter.geoms) > 0):
								flag = 1
						elif(isinstance(inter,Point)):
							if(len(inter.coords) > 0):
								flag = 1
						elif(isinstance(inter,MultiLineString)):
							if(len(inter.geoms) > 0):
								flag = 1
						if(flag == 0):
							point_height = interpolate_two_skel(isopter_labels_dict, distance_iso, distance_skel, iso, below_iso )
							break
					if(point_height == 0):
					        point_height = isopter_labels_dict[below_iso]['skeleton'][2][0]


				## else, if above_iso is specified, two point between the below and above (i.e. len(children) == 1)
				elif (above_iso != ''):
					if (len(isopter_children) > 0):
						#print("CASE2A")
						distance_below, point_below = find_distance_to_isopter(isopter_labels_dict, circle_point, below_iso)
						distance_above, point_above = find_distance_to_isopter(isopter_labels_dict, circle_point, above_iso)
						point_height = interpolate_two(isopter_labels_dict, distance_above, distance_below, above_iso, below_iso)
					else:
						## Isopter with only a scotoma child within it
						#print("CASE2B")
						distance_below, point_below = find_distance_to_isopter(isopter_labels_dict, circle_point, below_iso)
						distance_above, point_above = find_distance_to_isopter(isopter_labels_dict, circle_point, above_iso)
						min_distance_skel, min_point_skel = find_distance_to_skeleton(isopter_labels_dict, circle_point, below_iso)
						min_iso = below_iso
						min_distance = distance_below
						if(distance_above < distance_below):
							min_iso = above_iso
							min_distance = distance_above
						point_height =  interpolate_two_skel(isopter_labels_dict, min_distance, min_distance_skel, min_iso, below_iso)

				## else, skeleton either highest iso, or lowest scotoma
				else:
					#print("CASE3")
					#print("DEBUG: no above_iso, below_iso = %s" % below_iso)
					min_distance_below, min_point_below = find_distance_to_isopter(isopter_labels_dict, circle_point, below_iso)
					min_distance_skel, min_point_skel = find_distance_to_skeleton(isopter_labels_dict, circle_point, below_iso)
					point_height =  interpolate_two_skel(isopter_labels_dict, min_distance_below, min_distance_skel, below_iso, below_iso)

			## Flip the Y
			#my_point_temp = circle_point[1] - mY
			#my_point_temp = mY - my_point_temp
			#test_point = Point(circle_point[0], my_point_temp)
		

			## Add point and height to exiting lists
			points.append(test_point)
			heights.append(point_height * z_scale)
			volume = area*point_height
			volumes.append(volume)
			total_volume = total_volume + volume

			#print("DEBUG: r=%d, theta=%d, X=%0.2f, Y=%0.2f, z=%0.2f" % (i,j,test_point.x, test_point.y, point_height))


			## The rest of the code at in this scope should be removed.
			## Except for the area loop at the end.
			## And the calculation of the volume will need to re-integrated



			
			## Determine areas based upon point-height, and accumulate into isopter_area dictionary
			for iso in isopter_height.keys():
				if point_height > isopter_height[iso]: 
					isopter_area[iso] += area



	## Print out volume -- eventually this should be printed out to the METADATA/INFO file
	metadata_file = re.sub(".html$", ".summary", output_file)
	#print ('Metadata file is ', metadata_file)
	f = open(metadata_file, "w")
	f.write("Volume\t%f\n" % total_volume)
	f.write("Laterality\t%s\n" % eye)
	## Now print the areas of each isopter to the metadata_file
	for iso in isopter_height.keys():
		f.write("%s\t%s\n" % (iso, isopter_area[iso]))
	f.close()

	#print (total_volume)

	## DEBUG
	#print(heights)

	if(len(points) != len(heights)):
		print("# points (%d) != # heights (%d)" % (len(points), len(heights)) )

	## Generate the 3D graph
	#graph(circle_points_list, heights, isopter_labels_dict, output_file, z_scale)
	graph(points, heights, isopter_labels_dict, intersects_dict, output_file, z_scale, mX, mY, tX, tY, bX, bY, rX, rY, lX, lY, degree_len)


def skeleton_z(isopter_labels_dict, z_scale, each):
	#print (isopter_labels_dict[each]['name'])
	z_height = 0
	if "I1e" in isopter_labels_dict[each]['name']:
		z_height = 4.25
		#z_height = 127.5
		if "sc" in isopter_labels_dict[each]['name']:
			z_height = 3.75
			#z_height = 112.5
	if "I2e" in isopter_labels_dict[each]['name']:
		z_height = 3.75
		#z_height = 112.5
		if "sc" in isopter_labels_dict[each]['name']:
			z_height = 3.25
			#z_height = 97.5
	if "I3e" in isopter_labels_dict[each]['name']:
		z_height = 3.25
		#z_height = 97.5
		if "sc" in isopter_labels_dict[each]['name']:
			z_height = 2.75
			#z_height = 82.5
	if "I4e" in isopter_labels_dict[each]['name']:
		z_height = 2.75
		#z_height = 82.5
		if "sc" in isopter_labels_dict[each]['name']:
			z_height = 2.0
			#z_height = 60
	if "III4e" in isopter_labels_dict[each]['name']:
		z_height = 2.0
		#z_height = 60
		if "sc" in isopter_labels_dict[each]['name']:
			z_height = 1.0
			#z_height = 30
	if "V4e" in isopter_labels_dict[each]['name']:
		z_height = 1.0
		#z_height = 30
		if "sc" in isopter_labels_dict[each]['name']:
			z_height = 0.25
			#z_height = 7.5

	return (z_height)

def make_skeletons(isopter_labels_dict, each, scotoma_children, isopter_children, lX, rX, tY, bY):
	#Need to go through and find if no children, or only scotoma children
	#Going to test with I2e sc

	max_x, max_y = find_max_x_y(isopter_labels_dict, rX, bY)

#	print ("max_x", max_x)
#	print ("max_y", max_y)
	
	max_x = int(max_x)
	max_y = int(max_y)

	a = np.zeros((max_x, max_y)).astype(int)

	for i in range(max_x):
		for j in range(max_y):
			matrix_Point = Point(i,j)
			point_test  = isopter_labels_dict[each]['polygon'].contains(matrix_Point)
			if point_test:
				a[i][j] = 1
	if len(scotoma_children) > 0:
		for child in scotoma_children:
			for i in range(max_x):
				for j in range(max_y):
					matrix_Point = Point(i,j)
					point_test  = isopter_labels_dict[child]['polygon'].contains(matrix_Point)
					if point_test:
						a[i][j] = 0
	if len(isopter_children) > 0:
		for child in isopter_children:
			for i in range(max_x):
				for j in range(max_y):
					matrix_Point = Point(i,j)
					point_test  = isopter_labels_dict[child]['polygon'].contains(matrix_Point)
					if point_test:
						a[i][j] = 0
	skeleton = skeletonize(a)
	skeleton_x = []
	skeleton_y = []
	for i in range(max_x):
		for j in range(max_y):
			if skeleton[i][j] == 1:
				skeleton_x.append(i)
				skeleton_y.append(j)

	return (skeleton_x, skeleton_y, a, skeleton)

def display_skeleton(a, skeleton):

	# display results
	fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), sharex=True, sharey=True)
	ax = axes.ravel()

	ax[0].imshow(a, cmap=plt.cm.gray)
	ax[0].axis('off')
	ax[0].set_title('int_A', fontsize=20)		

	ax[1].imshow(skeleton, cmap=plt.cm.gray)
	ax[1].axis('off')
	ax[1].set_title('skeleton', fontsize=20)
	fig.tight_layout()

	plt.show()


def find_max_x_y(isopter_labels_dict, rX, bY):
	max_x = 0
	max_y = 0
	x_points = []
	y_points = []
	for each in isopter_labels_dict:
		points = isopter_labels_dict[each]['polygon'].exterior.coords
		for point in points:
			x_points.append(point[0])
			y_points.append(point[1])
	max_x = max(x_points)
	max_y = max(y_points)
	if max_x < rX:
		max_x = rX
	if max_y < bY:
		max_y = bY
	max_x = max_x + 50
	max_y = max_y + 50
	return (max_x, max_y)	

def graph(points_list, heights_list, isopter_labels_dict, intersects_dict, output_file, z_scale, mX, mY, tX, tY, bX, bY, rX, rY, lX, lY, degree_len):

	x_points_spline = []
	y_points_spline = []
	z_points_spline = []
	for each in isopter_labels_dict:
		points = isopter_labels_dict[each]['polygon'].exterior.coords
		#print(isopter_labels_dict[each]['name'])
		for point in points:
			x_points_spline.append(point[0])
			y_points_spline.append(point[1])
			z_points_spline.append(z_scale * isopter_labels_dict[each]['height'])
		## Only add the skeleton points IF there are no children of each
		if ('skeleton' in isopter_labels_dict[each]):
			iso_children, sc_children = find_children(isopter_labels_dict, intersects_dict, each)
			all_children = iso_children + sc_children
			if len(all_children) == 0:
				num_skel_points = len(isopter_labels_dict[each]['skeleton'][0])
				for i in range(0,num_skel_points):
					#print (isopter_labels_dict[each]['skeleton'][0][i], isopter_labels_dict[each]['skeleton'][1][i])
					x_points_spline.append(isopter_labels_dict[each]['skeleton'][0][i])
					y_points_spline.append(isopter_labels_dict[each]['skeleton'][1][i])
					z_points_spline.append(z_scale * isopter_labels_dict[each]['skeleton'][2][i])

	for ii in range(0,len(points_list),1):
		x_points_spline.append(points_list[ii].x)
		y_points_spline.append(points_list[ii].y)
		z_points_spline.append(heights_list[ii])

	
	#print (go.Mesh3d(x=x_points_spline, y=y_points_spline, z=z_points_spline))
	
	## Figure out how to specify the RGB/RGBA colors!
	## plotly.express.colors.sequential.Viridis ==>
	## ['#440154', '#482878', '#3e4989', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725']

	my_viridis = [ [0.00000, 'rgba(0.00000,0.00000,0.00000,0.0)' ], [0.01000, 'rgba(0.26667,0.00392,0.32941,1.0)' ], [0.11111, 'rgba(0.28235,0.15686,0.47059,1.0)' ], [0.22222, 'rgba(0.24314,0.28627,0.53725,1.0)' ], [0.33333, 'rgba(0.19216,0.40784,0.55686,1.0)' ], [0.44444, 'rgba(0.14902,0.50980,0.55686,1.0)' ], [0.55555, 'rgba(0.12157,0.61961,0.53725,1.0)' ], [0.66666, 'rgba(0.20784,0.71765,0.47451,1.0)' ], [0.77777, 'rgba(0.43137,0.80784,0.34510,1.0)' ], [0.88888, 'rgba(0.70980,0.87059,0.16863,1.0)' ], [1.00000, 'rgba(0.99216,0.90588,0.14510,1.0)' ] ]
	## Build a TURBO equivalent, with a transparent lowest color

	my_turbo = [ [0.00000,'rgba(0.18995,0.07176,0.23217,0.0)'], [ 0.002, 'rgba(0.19483,0.08339,0.26149,1.0)'], [ 0.00784, 'rgba(0.19956,0.09498,0.29024,1.0)'], [ 0.01176, 'rgba(0.20415,0.10652,0.31844,1.0)'], [ 0.01569, 'rgba(0.20860,0.11802,0.34607,1.0)'], [ 0.01961, 'rgba(0.21291,0.12947,0.37314,1.0)'], [ 0.02353, 'rgba(0.21708,0.14087,0.39964,1.0)'], [ 0.02745, 'rgba(0.22111,0.15223,0.42558,1.0)'], [ 0.03137, 'rgba(0.22500,0.16354,0.45096,1.0)'], [ 0.03529, 'rgba(0.22875,0.17481,0.47578,1.0)'], [ 0.03922, 'rgba(0.23236,0.18603,0.50004,1.0)'], [ 0.04314, 'rgba(0.23582,0.19720,0.52373,1.0)'], [ 0.04706, 'rgba(0.23915,0.20833,0.54686,1.0)'], [ 0.05098, 'rgba(0.24234,0.21941,0.56942,1.0)'], [ 0.05490, 'rgba(0.24539,0.23044,0.59142,1.0)'], [ 0.05882, 'rgba(0.24830,0.24143,0.61286,1.0)'], [ 0.06275, 'rgba(0.25107,0.25237,0.63374,1.0)'], [ 0.06667, 'rgba(0.25369,0.26327,0.65406,1.0)'], [ 0.07059, 'rgba(0.25618,0.27412,0.67381,1.0)'], [ 0.07451, 'rgba(0.25853,0.28492,0.69300,1.0)'], [ 0.07843, 'rgba(0.26074,0.29568,0.71162,1.0)'], [ 0.08235, 'rgba(0.26280,0.30639,0.72968,1.0)'], [ 0.08627, 'rgba(0.26473,0.31706,0.74718,1.0)'], [ 0.09020, 'rgba(0.26652,0.32768,0.76412,1.0)'], [ 0.09412, 'rgba(0.26816,0.33825,0.78050,1.0)'], [ 0.09804, 'rgba(0.26967,0.34878,0.79631,1.0)'], [ 0.10196, 'rgba(0.27103,0.35926,0.81156,1.0)'], [ 0.10588, 'rgba(0.27226,0.36970,0.82624,1.0)'], [ 0.10980, 'rgba(0.27334,0.38008,0.84037,1.0)'], [ 0.11373, 'rgba(0.27429,0.39043,0.85393,1.0)'], [ 0.11765, 'rgba(0.27509,0.40072,0.86692,1.0)'], [ 0.12157, 'rgba(0.27576,0.41097,0.87936,1.0)'], [ 0.12549, 'rgba(0.27628,0.42118,0.89123,1.0)'], [ 0.12941, 'rgba(0.27667,0.43134,0.90254,1.0)'], [ 0.13333, 'rgba(0.27691,0.44145,0.91328,1.0)'], [ 0.13725, 'rgba(0.27701,0.45152,0.92347,1.0)'], [ 0.14118, 'rgba(0.27698,0.46153,0.93309,1.0)'], [ 0.14510, 'rgba(0.27680,0.47151,0.94214,1.0)'], [ 0.14902, 'rgba(0.27648,0.48144,0.95064,1.0)'], [ 0.15294, 'rgba(0.27603,0.49132,0.95857,1.0)'], [ 0.15686, 'rgba(0.27543,0.50115,0.96594,1.0)'], [ 0.16078, 'rgba(0.27469,0.51094,0.97275,1.0)'], [ 0.16471, 'rgba(0.27381,0.52069,0.97899,1.0)'], [ 0.16863, 'rgba(0.27273,0.53040,0.98461,1.0)'], [ 0.17255, 'rgba(0.27106,0.54015,0.98930,1.0)'], [ 0.17647, 'rgba(0.26878,0.54995,0.99303,1.0)'], [ 0.18039, 'rgba(0.26592,0.55979,0.99583,1.0)'], [ 0.18431, 'rgba(0.26252,0.56967,0.99773,1.0)'], [ 0.18824, 'rgba(0.25862,0.57958,0.99876,1.0)'], [ 0.19216, 'rgba(0.25425,0.58950,0.99896,1.0)'], [ 0.19608, 'rgba(0.24946,0.59943,0.99835,1.0)'], [ 0.20000, 'rgba(0.24427,0.60937,0.99697,1.0)'], [ 0.20392, 'rgba(0.23874,0.61931,0.99485,1.0)'], [ 0.20784, 'rgba(0.23288,0.62923,0.99202,1.0)'], [ 0.21176, 'rgba(0.22676,0.63913,0.98851,1.0)'], [ 0.21569, 'rgba(0.22039,0.64901,0.98436,1.0)'], [ 0.21961, 'rgba(0.21382,0.65886,0.97959,1.0)'], [ 0.22353, 'rgba(0.20708,0.66866,0.97423,1.0)'], [ 0.22745, 'rgba(0.20021,0.67842,0.96833,1.0)'], [ 0.23137, 'rgba(0.19326,0.68812,0.96190,1.0)'], [ 0.23529, 'rgba(0.18625,0.69775,0.95498,1.0)'], [ 0.23922, 'rgba(0.17923,0.70732,0.94761,1.0)'], [ 0.24314, 'rgba(0.17223,0.71680,0.93981,1.0)'], [ 0.24706, 'rgba(0.16529,0.72620,0.93161,1.0)'], [ 0.25098, 'rgba(0.15844,0.73551,0.92305,1.0)'], [ 0.25490, 'rgba(0.15173,0.74472,0.91416,1.0)'], [ 0.25882, 'rgba(0.14519,0.75381,0.90496,1.0)'], [ 0.26275, 'rgba(0.13886,0.76279,0.89550,1.0)'], [ 0.26667, 'rgba(0.13278,0.77165,0.88580,1.0)'], [ 0.27059, 'rgba(0.12698,0.78037,0.87590,1.0)'], [ 0.27451, 'rgba(0.12151,0.78896,0.86581,1.0)'], [ 0.27843, 'rgba(0.11639,0.79740,0.85559,1.0)'], [ 0.28235, 'rgba(0.11167,0.80569,0.84525,1.0)'], [ 0.28627, 'rgba(0.10738,0.81381,0.83484,1.0)'], [ 0.29020, 'rgba(0.10357,0.82177,0.82437,1.0)'], [ 0.29412, 'rgba(0.10026,0.82955,0.81389,1.0)'], [ 0.29804, 'rgba(0.09750,0.83714,0.80342,1.0)'], [ 0.30196, 'rgba(0.09532,0.84455,0.79299,1.0)'], [ 0.30588, 'rgba(0.09377,0.85175,0.78264,1.0)'], [ 0.30980, 'rgba(0.09287,0.85875,0.77240,1.0)'], [ 0.31373, 'rgba(0.09267,0.86554,0.76230,1.0)'], [ 0.31765, 'rgba(0.09320,0.87211,0.75237,1.0)'], [ 0.32157, 'rgba(0.09451,0.87844,0.74265,1.0)'], [ 0.32549, 'rgba(0.09662,0.88454,0.73316,1.0)'], [ 0.32941, 'rgba(0.09958,0.89040,0.72393,1.0)'], [ 0.33333, 'rgba(0.10342,0.89600,0.71500,1.0)'], [ 0.33725, 'rgba(0.10815,0.90142,0.70599,1.0)'], [ 0.34118, 'rgba(0.11374,0.90673,0.69651,1.0)'], [ 0.34510, 'rgba(0.12014,0.91193,0.68660,1.0)'], [ 0.34902, 'rgba(0.12733,0.91701,0.67627,1.0)'], [ 0.35294, 'rgba(0.13526,0.92197,0.66556,1.0)'], [ 0.35686, 'rgba(0.14391,0.92680,0.65448,1.0)'], [ 0.36078, 'rgba(0.15323,0.93151,0.64308,1.0)'], [ 0.36471, 'rgba(0.16319,0.93609,0.63137,1.0)'], [ 0.36863, 'rgba(0.17377,0.94053,0.61938,1.0)'], [ 0.37255, 'rgba(0.18491,0.94484,0.60713,1.0)'], [ 0.37647, 'rgba(0.19659,0.94901,0.59466,1.0)'], [ 0.38039, 'rgba(0.20877,0.95304,0.58199,1.0)'], [ 0.38431, 'rgba(0.22142,0.95692,0.56914,1.0)'], [ 0.38824, 'rgba(0.23449,0.96065,0.55614,1.0)'], [ 0.39216, 'rgba(0.24797,0.96423,0.54303,1.0)'], [ 0.39608, 'rgba(0.26180,0.96765,0.52981,1.0)'], [ 0.40000, 'rgba(0.27597,0.97092,0.51653,1.0)'], [ 0.40392, 'rgba(0.29042,0.97403,0.50321,1.0)'], [ 0.40784, 'rgba(0.30513,0.97697,0.48987,1.0)'], [ 0.41176, 'rgba(0.32006,0.97974,0.47654,1.0)'], [ 0.41569, 'rgba(0.33517,0.98234,0.46325,1.0)'], [ 0.41961, 'rgba(0.35043,0.98477,0.45002,1.0)'], [ 0.42353, 'rgba(0.36581,0.98702,0.43688,1.0)'], [ 0.42745, 'rgba(0.38127,0.98909,0.42386,1.0)'], [ 0.43137, 'rgba(0.39678,0.99098,0.41098,1.0)'], [ 0.43529, 'rgba(0.41229,0.99268,0.39826,1.0)'], [ 0.43922, 'rgba(0.42778,0.99419,0.38575,1.0)'], [ 0.44314, 'rgba(0.44321,0.99551,0.37345,1.0)'], [ 0.44706, 'rgba(0.45854,0.99663,0.36140,1.0)'], [ 0.45098, 'rgba(0.47375,0.99755,0.34963,1.0)'], [ 0.45490, 'rgba(0.48879,0.99828,0.33816,1.0)'], [ 0.45882, 'rgba(0.50362,0.99879,0.32701,1.0)'], [ 0.46275, 'rgba(0.51822,0.99910,0.31622,1.0)'], [ 0.46667, 'rgba(0.53255,0.99919,0.30581,1.0)'], [ 0.47059, 'rgba(0.54658,0.99907,0.29581,1.0)'], [ 0.47451, 'rgba(0.56026,0.99873,0.28623,1.0)'], [ 0.47843, 'rgba(0.57357,0.99817,0.27712,1.0)'], [ 0.48235, 'rgba(0.58646,0.99739,0.26849,1.0)'], [ 0.48627, 'rgba(0.59891,0.99638,0.26038,1.0)'], [ 0.49020, 'rgba(0.61088,0.99514,0.25280,1.0)'], [ 0.49412, 'rgba(0.62233,0.99366,0.24579,1.0)'], [ 0.49804, 'rgba(0.63323,0.99195,0.23937,1.0)'], [ 0.50196, 'rgba(0.64362,0.98999,0.23356,1.0)'], [ 0.50588, 'rgba(0.65394,0.98775,0.22835,1.0)'], [ 0.50980, 'rgba(0.66428,0.98524,0.22370,1.0)'], [ 0.51373, 'rgba(0.67462,0.98246,0.21960,1.0)'], [ 0.51765, 'rgba(0.68494,0.97941,0.21602,1.0)'], [ 0.52157, 'rgba(0.69525,0.97610,0.21294,1.0)'], [ 0.52549, 'rgba(0.70553,0.97255,0.21032,1.0)'], [ 0.52941, 'rgba(0.71577,0.96875,0.20815,1.0)'], [ 0.53333, 'rgba(0.72596,0.96470,0.20640,1.0)'], [ 0.53725, 'rgba(0.73610,0.96043,0.20504,1.0)'], [ 0.54118, 'rgba(0.74617,0.95593,0.20406,1.0)'], [ 0.54510, 'rgba(0.75617,0.95121,0.20343,1.0)'], [ 0.54902, 'rgba(0.76608,0.94627,0.20311,1.0)'], [ 0.55294, 'rgba(0.77591,0.94113,0.20310,1.0)'], [ 0.55686, 'rgba(0.78563,0.93579,0.20336,1.0)'], [ 0.56078, 'rgba(0.79524,0.93025,0.20386,1.0)'], [ 0.56471, 'rgba(0.80473,0.92452,0.20459,1.0)'], [ 0.56863, 'rgba(0.81410,0.91861,0.20552,1.0)'], [ 0.57255, 'rgba(0.82333,0.91253,0.20663,1.0)'], [ 0.57647, 'rgba(0.83241,0.90627,0.20788,1.0)'], [ 0.58039, 'rgba(0.84133,0.89986,0.20926,1.0)'], [ 0.58431, 'rgba(0.85010,0.89328,0.21074,1.0)'], [ 0.58824, 'rgba(0.85868,0.88655,0.21230,1.0)'], [ 0.59216, 'rgba(0.86709,0.87968,0.21391,1.0)'], [ 0.59608, 'rgba(0.87530,0.87267,0.21555,1.0)'], [ 0.60000, 'rgba(0.88331,0.86553,0.21719,1.0)'], [ 0.60392, 'rgba(0.89112,0.85826,0.21880,1.0)'], [ 0.60784, 'rgba(0.89870,0.85087,0.22038,1.0)'], [ 0.61176, 'rgba(0.90605,0.84337,0.22188,1.0)'], [ 0.61569, 'rgba(0.91317,0.83576,0.22328,1.0)'], [ 0.61961, 'rgba(0.92004,0.82806,0.22456,1.0)'], [ 0.62353, 'rgba(0.92666,0.82025,0.22570,1.0)'], [ 0.62745, 'rgba(0.93301,0.81236,0.22667,1.0)'], [ 0.63137, 'rgba(0.93909,0.80439,0.22744,1.0)'], [ 0.63529, 'rgba(0.94489,0.79634,0.22800,1.0)'], [ 0.63922, 'rgba(0.95039,0.78823,0.22831,1.0)'], [ 0.64314, 'rgba(0.95560,0.78005,0.22836,1.0)'], [ 0.64706, 'rgba(0.96049,0.77181,0.22811,1.0)'], [ 0.65098, 'rgba(0.96507,0.76352,0.22754,1.0)'], [ 0.65490, 'rgba(0.96931,0.75519,0.22663,1.0)'], [ 0.65882, 'rgba(0.97323,0.74682,0.22536,1.0)'], [ 0.66275, 'rgba(0.97679,0.73842,0.22369,1.0)'], [ 0.66667, 'rgba(0.98000,0.73000,0.22161,1.0)'], [ 0.67059, 'rgba(0.98289,0.72140,0.21918,1.0)'], [ 0.67451, 'rgba(0.98549,0.71250,0.21650,1.0)'], [ 0.67843, 'rgba(0.98781,0.70330,0.21358,1.0)'], [ 0.68235, 'rgba(0.98986,0.69382,0.21043,1.0)'], [ 0.68627, 'rgba(0.99163,0.68408,0.20706,1.0)'], [ 0.69020, 'rgba(0.99314,0.67408,0.20348,1.0)'], [ 0.69412, 'rgba(0.99438,0.66386,0.19971,1.0)'], [ 0.69804, 'rgba(0.99535,0.65341,0.19577,1.0)'], [ 0.70196, 'rgba(0.99607,0.64277,0.19165,1.0)'], [ 0.70588, 'rgba(0.99654,0.63193,0.18738,1.0)'], [ 0.70980, 'rgba(0.99675,0.62093,0.18297,1.0)'], [ 0.71373, 'rgba(0.99672,0.60977,0.17842,1.0)'], [ 0.71765, 'rgba(0.99644,0.59846,0.17376,1.0)'], [ 0.72157, 'rgba(0.99593,0.58703,0.16899,1.0)'], [ 0.72549, 'rgba(0.99517,0.57549,0.16412,1.0)'], [ 0.72941, 'rgba(0.99419,0.56386,0.15918,1.0)'], [ 0.73333, 'rgba(0.99297,0.55214,0.15417,1.0)'], [ 0.73725, 'rgba(0.99153,0.54036,0.14910,1.0)'], [ 0.74118, 'rgba(0.98987,0.52854,0.14398,1.0)'], [ 0.74510, 'rgba(0.98799,0.51667,0.13883,1.0)'], [ 0.74902, 'rgba(0.98590,0.50479,0.13367,1.0)'], [ 0.75294, 'rgba(0.98360,0.49291,0.12849,1.0)'], [ 0.75686, 'rgba(0.98108,0.48104,0.12332,1.0)'], [ 0.76078, 'rgba(0.97837,0.46920,0.11817,1.0)'], [ 0.76471, 'rgba(0.97545,0.45740,0.11305,1.0)'], [ 0.76863, 'rgba(0.97234,0.44565,0.10797,1.0)'], [ 0.77255, 'rgba(0.96904,0.43399,0.10294,1.0)'], [ 0.77647, 'rgba(0.96555,0.42241,0.09798,1.0)'], [ 0.78039, 'rgba(0.96187,0.41093,0.09310,1.0)'], [ 0.78431, 'rgba(0.95801,0.39958,0.08831,1.0)'], [ 0.78824, 'rgba(0.95398,0.38836,0.08362,1.0)'], [ 0.79216, 'rgba(0.94977,0.37729,0.07905,1.0)'], [ 0.79608, 'rgba(0.94538,0.36638,0.07461,1.0)'], [ 0.80000, 'rgba(0.94084,0.35566,0.07031,1.0)'], [ 0.80392, 'rgba(0.93612,0.34513,0.06616,1.0)'], [ 0.80784, 'rgba(0.93125,0.33482,0.06218,1.0)'], [ 0.81176, 'rgba(0.92623,0.32473,0.05837,1.0)'], [ 0.81569, 'rgba(0.92105,0.31489,0.05475,1.0)'], [ 0.81961, 'rgba(0.91572,0.30530,0.05134,1.0)'], [ 0.82353, 'rgba(0.91024,0.29599,0.04814,1.0)'], [ 0.82745, 'rgba(0.90463,0.28696,0.04516,1.0)'], [ 0.83137, 'rgba(0.89888,0.27824,0.04243,1.0)'], [ 0.83529, 'rgba(0.89298,0.26981,0.03993,1.0)'], [ 0.83922, 'rgba(0.88691,0.26152,0.03753,1.0)'], [ 0.84314, 'rgba(0.88066,0.25334,0.03521,1.0)'], [ 0.84706, 'rgba(0.87422,0.24526,0.03297,1.0)'], [ 0.85098, 'rgba(0.86760,0.23730,0.03082,1.0)'], [ 0.85490, 'rgba(0.86079,0.22945,0.02875,1.0)'], [ 0.85882, 'rgba(0.85380,0.22170,0.02677,1.0)'], [ 0.86275, 'rgba(0.84662,0.21407,0.02487,1.0)'], [ 0.86667, 'rgba(0.83926,0.20654,0.02305,1.0)'], [ 0.87059, 'rgba(0.83172,0.19912,0.02131,1.0)'], [ 0.87451, 'rgba(0.82399,0.19182,0.01966,1.0)'], [ 0.87843, 'rgba(0.81608,0.18462,0.01809,1.0)'], [ 0.88235, 'rgba(0.80799,0.17753,0.01660,1.0)'], [ 0.88627, 'rgba(0.79971,0.17055,0.01520,1.0)'], [ 0.89020, 'rgba(0.79125,0.16368,0.01387,1.0)'], [ 0.89412, 'rgba(0.78260,0.15693,0.01264,1.0)'], [ 0.89804, 'rgba(0.77377,0.15028,0.01148,1.0)'], [ 0.90196, 'rgba(0.76476,0.14374,0.01041,1.0)'], [ 0.90588, 'rgba(0.75556,0.13731,0.00942,1.0)'], [ 0.90980, 'rgba(0.74617,0.13098,0.00851,1.0)'], [ 0.91373, 'rgba(0.73661,0.12477,0.00769,1.0)'], [ 0.91765, 'rgba(0.72686,0.11867,0.00695,1.0)'], [ 0.92157, 'rgba(0.71692,0.11268,0.00629,1.0)'], [ 0.92549, 'rgba(0.70680,0.10680,0.00571,1.0)'], [ 0.92941, 'rgba(0.69650,0.10102,0.00522,1.0)'], [ 0.93333, 'rgba(0.68602,0.09536,0.00481,1.0)'], [ 0.93725, 'rgba(0.67535,0.08980,0.00449,1.0)'], [ 0.94118, 'rgba(0.66449,0.08436,0.00424,1.0)'], [ 0.94510, 'rgba(0.65345,0.07902,0.00408,1.0)'], [ 0.94902, 'rgba(0.64223,0.07380,0.00401,1.0)'], [ 0.95294, 'rgba(0.63082,0.06868,0.00401,1.0)'], [ 0.95686, 'rgba(0.61923,0.06367,0.00410,1.0)'], [ 0.96078, 'rgba(0.60746,0.05878,0.00427,1.0)'], [ 0.96471, 'rgba(0.59550,0.05399,0.00453,1.0)'], [ 0.96863, 'rgba(0.58336,0.04931,0.00486,1.0)'], [ 0.97255, 'rgba(0.57103,0.04474,0.00529,1.0)'], [ 0.97647, 'rgba(0.55852,0.04028,0.00579,1.0)'], [ 0.98039, 'rgba(0.54583,0.03593,0.00638,1.0)'], [ 0.98431, 'rgba(0.53295,0.03169,0.00705,1.0)'], [ 0.98824, 'rgba(0.51989,0.02756,0.00780,1.0)'], [ 0.99216, 'rgba(0.50664,0.02354,0.00863,1.0)'], [ 0.99608, 'rgba(0.49321,0.01963,0.00955,1.0)'], [ 1.00000, 'rgba(0.47960,0.01583,0.01055) '] ]


	#fig = go.Figure(data=[go.Mesh3d(x=x_points_spline, y=y_points_spline, z=z_points_spline, colorscale = my_turbo , cmin=0.5*z_scale, cmax=4.5*z_scale, intensity=z_points_spline, opacity=0.50)])
	#fig = go.Figure(data=[go.Mesh3d(x=x_points_spline, y=y_points_spline, z=z_points_spline, colorscale = my_viridis , cmin=0.0*z_scale, cmax=4.5*z_scale, intensity=z_points_spline, opacity=0.50)])
	fig = go.Figure(data=[go.Mesh3d(x=x_points_spline, y=y_points_spline, z=z_points_spline, colorscale = my_viridis , cmin=0.0*z_scale, cmax=4.5*z_scale, intensity=z_points_spline, opacity=0.90, flatshading=False)])
	#fig = go.Figure(data=[go.Mesh3d(x=x_points_spline, y=y_points_spline, z=z_points_spline, colorscale = "Viridis", cmin=0.5*z_scale, cmax=4.5*z_scale, intensity=z_points_spline, opacity=0.50)])

	## Need to flip the Y axis, as the original scoring is from the iPad/TruthMarker, with (0,0) in the upper left-hand corner.
	#fig.update_scenes(yaxis_autorange="reversed")
	#fig.update_scenes(projection_type="orthographic")
	fig.layout.scene.camera.projection.type = "orthographic"
	fig.layout.showlegend = False

	fig.update_layout(scene = dict(
		xaxis = dict(
			showbackground = False,
			showaxeslabels = False,
			visible = False,
			showticklabels = False),
		yaxis = dict(
			showbackground = False,
			showaxeslabels = False,
			visible = False,
			showticklabels = False),
		zaxis = dict(
			showbackground = False,
			showticklabels = False,
			visible = False,
			showaxeslabels = False)))


#	fig.add_trace(go.Scatter3d(x=x_points,y=y_points,z=z_points, mode="markers",marker=dict(color="grey")))

	x_points = { "V4e": [], "III4e": [], "I4e": [], "I3e": [], "I2e": [], "I1e": [] }
	y_points = { "V4e": [], "III4e": [], "I4e": [], "I3e": [], "I2e": [], "I1e": [] }
	z_points = { "V4e": [], "III4e": [], "I4e": [], "I3e": [], "I2e": [], "I1e": [] }
	for iso_label in isopter_labels_dict.keys():
		points = isopter_labels_dict[iso_label]['polygon'].exterior.coords
		type = isopter_labels_dict[iso_label]['name']
		height = 0
		if(type.startswith("V4e")):
			type = "V4e"
			height = 0.5 * z_scale
		elif(type.startswith("III4e")):
			type = "III4e"
			height = 1.5 * z_scale
		elif(type.startswith("I4e")):
			type = "I4e"
			height = 2.5 * z_scale
		elif(type.startswith("I3e")):
			type = "I3e"
			height = 3.0 * z_scale
		elif(type.startswith("I2e")):
			type = "I2e"
			height = 3.5 * z_scale
		elif(type.startswith("I1e")):
			type = "I1e"
			height = 4.0 * z_scale
		for point in points:
			x_points[type].append(point[0])
			y_points[type].append(point[1])
			z_points[type].append(height)
			
	fig.add_trace(go.Scatter3d(x=x_points["I1e"],y=y_points["I1e"],z=z_points["I1e"], mode="markers", marker=dict(color="green", size=2)))
	fig.add_trace(go.Scatter3d(x=x_points["I2e"],y=y_points["I2e"],z=z_points["I2e"], mode="markers", marker=dict(color="red", size=2)))
	fig.add_trace(go.Scatter3d(x=x_points["I3e"],y=y_points["I3e"],z=z_points["I3e"], mode="markers", marker=dict(color="gray", size=2)))
	fig.add_trace(go.Scatter3d(x=x_points["I4e"],y=y_points["I4e"],z=z_points["I4e"], mode="markers", marker=dict(color="blue", size=2)))
	fig.add_trace(go.Scatter3d(x=x_points["III4e"],y=y_points["III4e"],z=z_points["III4e"], mode="markers", marker=dict(color="brown", size=2)))
	fig.add_trace(go.Scatter3d(x=x_points["V4e"],y=y_points["V4e"],z=z_points["V4e"], mode="markers", marker=dict(color="purple", size=2)))
	## Add a point to compress the Z-scaling
	fig.add_trace(go.Scatter3d(x=[0],y=[0],z=[8], mode="markers", marker=dict(color='rgba(1,1,1,0)', size=2)))


	## Render the grid
	draw_grid(fig, degree_len, mX, mY)

	## Adding camera control - view from straight above (rather than default)
	camera = dict(eye=dict(x=0.0, y=0.0, z=1.0))
	fig.update_layout(scene_camera=camera)

	output_file_jpg = re.sub(".html$", ".jpg", output_file)
	fig.write_image(output_file_jpg)

	fig.write_html(output_file)
	#fig.add_trace(go.Scatter3d(x=[500],y=[500],z=[200]))

	#fig.update_layout(scene = dict(xaxis = dict(nticks=4, range=[0,800],),yaxis = dict(nticks=4, range=[0,800],),zaxis = dict(nticks=4, range=[-50,200],),))

def draw_grid(fig, degree_len, mX, mY):
	#Nicole is adding the following for the grid

	distance1 = math.tan(math.radians(15))*(70*degree_len)
	distance2 = math.tan(math.radians(30))*(70*degree_len)
	distance3 = math.tan(math.radians(165))*(74*degree_len)
	distance4 = math.tan(math.radians(150))*(74*degree_len)

	## Line at 15 degrees CW from vertical, from 5 degrees to 70 degrees
	x10_15 = (5*degree_len*(math.sin(math.radians(15)))) + mX
	y10_15 = (5*degree_len*(math.cos(math.radians(15)))) + mY
	y90_15 = (70*degree_len) + mY
	x90_15 = mX + distance1
	fig.add_trace(go.Scatter3d(x = [x10_15,x90_15], y = [y10_15,y90_15], z = [0,0], mode='lines', line=dict(color = 'black')))

	## Line at 30 degrees CW from vertical, from 5 degrees to 70 degrees
	x10_30 = (5*degree_len*(math.sin(math.radians(30)))) + mX
	y10_30 = (5*degree_len*(math.cos(math.radians(30)))) + mY
	x90_30 = mX + distance2
	y90_30 = (70*degree_len) + mY
	fig.add_trace(go.Scatter3d(x = [x10_30,x90_30], y = [y10_30,y90_30], z = [0,0], mode='lines', line=dict(color = 'black')))

	## Line at 150 degrees CW from vertical, from 5 degrees to 74 degrees
	x10_150 = (5*degree_len*(math.sin(math.radians(150)))) + mX
	y10_150 = (5*degree_len*(math.cos(math.radians(150)))) + mY
	y90_150 = mY - (74*degree_len)
	x90_150 = mX - distance4
	fig.add_trace(go.Scatter3d(x = [x10_150,x90_150], y = [y10_150,y90_150], z = [0,0], mode='lines', line=dict(color = 'black')))

	## Line at 165 degrees CW from vertical, from 5 degrees to 74 degrees
	x10_165 = (5*degree_len*(math.sin(math.radians(165)))) + mX
	y10_165 = (5*degree_len*(math.cos(math.radians(165)))) + mY
	y90_165 = mY - (74*degree_len)
	x90_165 = mX - distance3
	fig.add_trace(go.Scatter3d(x = [x10_165,x90_165], y = [y10_165,y90_165], z = [0,0], mode='lines', line=dict(color = 'black')))

	## Line at 195 degrees CW (165 degrees CCW) from vertical, from 5 degrees to 74 degrees
	x10_195 = (5*degree_len*(math.sin(math.radians(195)))) + mX
	y10_195 = (5*degree_len*(math.cos(math.radians(195)))) + mY
	y90_195 = mY - (74*degree_len)
	x90_195 = mX + distance3
	fig.add_trace(go.Scatter3d(x = [x10_195,x90_195], y = [y10_195,y90_195], z = [0,0], mode='lines', line=dict(color = 'black')))

	## Line at 210 degrees CW (150 degrees CCW) from vertical, from 5 degrees to 74 degrees
	x10_210 = (5*degree_len*(math.sin(math.radians(210)))) + mX
	y10_210 = (5*degree_len*(math.cos(math.radians(210)))) + mY
	y90_210 = mY - (74*degree_len)
	x90_210 = mX + distance4
	fig.add_trace(go.Scatter3d(x = [x10_210,x90_210], y = [y10_210,y90_210], z = [0,0], mode='lines', line=dict(color = 'black')))

	## Line at 330 degrees CW (30 degrees CCW) from vertical, from 5 degrees to 70 degrees
	x10_330 = (5*degree_len*(math.sin(math.radians(330)))) + mX
	y10_330 = (5*degree_len*(math.cos(math.radians(330)))) + mY
	x90_330 = mX - distance2
	y90_330 = (70*degree_len) + mY
	fig.add_trace(go.Scatter3d(x = [x10_330,x90_330], y = [y10_330,y90_330], z = [0,0], mode='lines', line=dict(color = 'black')))

	## Line at 345 degrees CW (15 degrees CCW) from vertical, from 5 degrees to 70 degrees
	x10_345 = (5*degree_len*(math.sin(math.radians(345)))) + mX
	y10_345 = (5*degree_len*(math.cos(math.radians(345)))) + mY
	x90_345 = mX - distance1
	y90_345 = (70*degree_len) + mY
	fig.add_trace(go.Scatter3d(x = [x10_345,x90_345], y = [y10_345,y90_345], z = [0,0], mode='lines', line=dict(color = 'black')))

	## Capstone line at 70 degrees at top -- 34.4 degrees "out"
	x_dist = 90 * degree_len * (math.sin(math.radians(38.95)))
	y1 = mY + (70*degree_len)
	x1 = mX + x_dist
	x2 = mX - x_dist
	fig.add_trace(go.Scatter3d(x = [x1,x2], y = [y1,y1], z = [0,0], mode='lines', line=dict(color = 'black', width = 5)))

	## Capstone line at 74 degrees at bottom -- 34.4 degrees "out"
	x_dist = 90 * degree_len * (math.sin(math.radians(34.7)))
	y1 = mY - (74*degree_len)
	x1 = mX + x_dist
	x2 = mX - x_dist
	fig.add_trace(go.Scatter3d(x = [x1,x2], y = [y1,y1], z = [0,0], mode='lines', line=dict(color = 'black', width = 5)))


	#x1 = (90*degree_len*(math.sin(math.radians(320.6)))) + mX
	#y1 = mY - (90*degree_len*(math.cos(math.radians(320.8))))
	#x2 = (90*degree_len*(math.sin(math.radians(39.2)))) + mX
	#y2 = mY - (90*degree_len*(math.cos(math.radians(39.2))))
	#fig.add_trace(go.Scatter3d(x = [x1,x2], y = [y1,y2], z = [0,0], mode='lines', line=dict(color = 'black', width = 5)))

	#x1 = (90*degree_len*(math.sin(math.radians(145.6)))) + mX
	#y = mY+(degree_len*74.3)
	#x2 = (90*degree_len*(math.sin(math.radians(214.4)))) + mX
	#fig.add_trace(go.Scatter3d(x = [x1,x2], y = [y,y], z = [0,0], mode='lines', line=dict(color = 'black', width = 5)))



	## All other radial lines from 5 drgrees to 90 degrees
	line_degree_list = [45,60,75,105,120,135,225,240,255,285,300,315]
	for each in line_degree_list:
		x5 = (5*degree_len*(math.sin(math.radians(each)))) + mX
		y5 = (5*degree_len*(math.cos(math.radians(each)))) + mY
		x90 = (90*degree_len*(math.sin(math.radians(each)))) + mX
		y90 = (90*degree_len*(math.cos(math.radians(each)))) + mY
		fig.add_trace(go.Scatter3d(x = [x5,x90], y = [y5,y90], z = [0,0], mode='lines', line=dict(color = 'black')))

	## Cardinal lines (horizontal and vertical)
	fig.add_trace(go.Scatter3d(x = [mX,mX], y = [(mY-(degree_len*74)),(mY+(degree_len*70))], z = [0,0], mode='lines', line=dict(color = 'black')))
	fig.add_trace(go.Scatter3d(x = [(mX+(90*degree_len)),(mX-(90*degree_len))], y = [mY,mY], z = [0,0], mode='lines', line=dict(color = 'black')))
	#fig.add_trace(go.Scatter3d(x = [mX,mX], y = [mY,(mY-(degree_len*74))], z = [0,0], mode='lines', line=dict(color = 'black')))
	#fig.add_trace(go.Scatter3d(x = [mX,mX], y = [mY,(mY+(degree_len*70))], z = [0,0], mode='lines', line=dict(color = 'black')))
	#fig.add_trace(go.Scatter3d(x = [mX,(mX-(90*degree_len))], y = [mY,mY], z = [0,0], mode='lines', line=dict(color = 'black')))
	#fig.add_trace(go.Scatter3d(x = [mX,(mX+(90*degree_len))], y = [mY,mY], z = [0,0], mode='lines', line=dict(color = 'black')))

	## Circles for distances from center with default thickness
	circle_deg_list = [5,10,20,40,50,70]
	for i in circle_deg_list:
		x_circle = []
		y_circle = []
		z_circle = []
		for j in range(0,361,1):
			x_circle.append(i*degree_len*(math.sin(math.radians(j)))+mX)
			y_circle.append(i*degree_len*(math.cos(math.radians(j)))+mY)
			z_circle.append(0)
		fig.add_trace(go.Scatter3d(x = x_circle, y = y_circle, z = z_circle, mode='lines', line=dict(color = 'black')))

	## Circles for distances from center with increased thickness
	circle_deg_list = [30,60]
	for i in circle_deg_list:
		x_circle = []
		y_circle = []
		z_circle = []
		for j in range(0,361,1):
			x_circle.append(i*degree_len*(math.sin(math.radians(j)))+mX)
			y_circle.append(i*degree_len*(math.cos(math.radians(j)))+mY)
			z_circle.append(0)
		fig.add_trace(go.Scatter3d(x = x_circle, y = y_circle, z = z_circle, mode='lines', line=dict(color = 'black', width = 5)))

	# 70 ==> 38.95
	# 74 ==> 34.7
	## Right-hand arcs for 80 and 90 degree lines
	x_circle = []
	y_circle = []
	z_circle = []
	for j in range(39,146,1):
		x_circle.append(90*degree_len*(math.sin(math.radians(j)))+mX)
		y_circle.append(90*degree_len*(math.cos(math.radians(j)))+mY)
		z_circle.append(0)
	x_circle.append(90*degree_len*(math.sin(math.radians(145.33)))+mX)
	y_circle.append(90*degree_len*(math.cos(math.radians(145.33)))+mY)
	z_circle.append(0)
	fig.add_trace(go.Scatter3d(x = x_circle, y = y_circle, z = z_circle, mode='lines', line=dict(color = 'black', width = 5)))


##
	x_circle = []
	y_circle = []
	z_circle = []
	for j in range(29,158,1):
		x_circle.append(80*degree_len*(math.sin(math.radians(j)))+mX)
		y_circle.append(80*degree_len*(math.cos(math.radians(j)))+mY)
		z_circle.append(0)
	x_circle.append(80*degree_len*(math.sin(math.radians(157.67)))+mX)
	y_circle.append(80*degree_len*(math.cos(math.radians(157.67)))+mY)
	z_circle.append(0)
	fig.add_trace(go.Scatter3d(x = x_circle, y = y_circle, z = z_circle, mode='lines', line=dict(color = 'black')))


	## Left-hand arcs for 80 and 90 degree lines
	x_circle = []
	y_circle = []
	z_circle = []
	x_circle.append(90*degree_len*(math.sin(math.radians(214.67)))+mX)
	y_circle.append(90*degree_len*(math.cos(math.radians(214.67)))+mY)
	z_circle.append(0)
	for j in range(215,322,1):
		x_circle.append(90*degree_len*(math.sin(math.radians(j)))+mX)
		y_circle.append(90*degree_len*(math.cos(math.radians(j)))+mY)
		z_circle.append(0)
	fig.add_trace(go.Scatter3d(x = x_circle, y = y_circle, z = z_circle, mode='lines', line=dict(color = 'black', width = 5)))
	x_circle = []
	y_circle = []
	z_circle = []
	x_circle.append(80*degree_len*(math.sin(math.radians(202.33)))+mX)
	y_circle.append(80*degree_len*(math.cos(math.radians(202.33)))+mY)
	z_circle.append(0)
	for j in range(203,332,1):
		x_circle.append(80*degree_len*(math.sin(math.radians(j)))+mX)
		y_circle.append(80*degree_len*(math.cos(math.radians(j)))+mY)
		z_circle.append(0)
	fig.add_trace(go.Scatter3d(x = x_circle, y = y_circle, z = z_circle, mode='lines', line=dict(color = 'black')))

	## 0 degrees is "up"
	## Add text labels
	fig.add_trace(go.Scatter3d(text = '90', x = [mX], y = [mY+((70+15)*degree_len*(math.cos(math.radians(0))))], z = [0], mode='text', textposition="middle center"))
	fig.add_trace(go.Scatter3d(text = '270', x = [mX], y = [mY+((70+15)*degree_len*(math.cos(math.radians(180))))], z = [0], mode='text', textposition="middle center"))
	fig.add_trace(go.Scatter3d(text = '0', x = [mX+((90+15)*degree_len*(math.sin(math.radians(90))))], y = [mY], z = [0], mode='text', textposition="middle center"))
	fig.add_trace(go.Scatter3d(text = '180', x = [mX+((90+15)*degree_len*(math.sin(math.radians(270))))], y = [mY], z = [0], mode='text', textposition="middle center"))




def interpolate_two(isopter_labels_dict, min_distance_above, min_distance_below, above_iso, below_iso):
	inverse_dist_1 = 1/min_distance_above
	inverse_dist_2 = 1/min_distance_below
	sum_inverse_dist = inverse_dist_1 + inverse_dist_2

	weight_1 = inverse_dist_1/sum_inverse_dist
	weight_2 = inverse_dist_2/sum_inverse_dist

	height_1 = isopter_labels_dict[above_iso]['height']
	height_2 = isopter_labels_dict[below_iso]['height']

#	print (height_1, weight_1)
#	print (height_2, weight_2)

	interpolated_height = height_1*weight_1 + height_2*weight_2

	return (interpolated_height)

def interpolate_two_skel(isopter_labels_dict, min_distance_non_skel, min_distance_skel, non_skel_iso, skel_iso):
	#print ("below_iso: ", below_iso)
	inverse_dist_1 = 1/min_distance_non_skel
	inverse_dist_2 = 1/min_distance_skel
	sum_inverse_dist = inverse_dist_1 + inverse_dist_2

	weight_1 = inverse_dist_1/sum_inverse_dist
	weight_2 = inverse_dist_2/sum_inverse_dist


	#print ("skeleton test")
	#print (isopter_labels_dict[below_iso]['skeleton'])

	height_1 = isopter_labels_dict[non_skel_iso]['height']
	height_2 = isopter_labels_dict[skel_iso]['skeleton'][2][0]

	#print (height_1*weight_1, height_2*weight_2)
	#print ((height_1*weight_1 + height_2*weight_2)/2)

	interpolated_height = (height_1*weight_1 + height_2*weight_2)

	return (interpolated_height)

#point_height = interpolate_three_skel(isopter_labels_dict, min_distance_c1, min_distance_c2, min_distance_skel, child_1, child_2, below_iso )
def interpolate_three_skel(isopter_labels_dict, min_distance_1, min_distance_2, min_distance_skel, non_skel_1, non_skel_2, skel_iso):
	#print ("below_iso: ", below_iso)
	inverse_dist_1 = 1/min_distance_1
	inverse_dist_2 = 1/min_distance_2
	inverse_dist_3 = 1/min_distance_skel
	sum_inverse_dist = inverse_dist_1 + inverse_dist_2 + inverse_dist_3

	weight_1 = inverse_dist_1/sum_inverse_dist
	weight_2 = inverse_dist_2/sum_inverse_dist
	weight_3 = inverse_dist_3/sum_inverse_dist

	#print ("skeleton test")
	#print (isopter_labels_dict[below_iso]['skeleton'])

	height_1 = isopter_labels_dict[non_skel_1]['height']
	height_2 = isopter_labels_dict[non_skel_2]['height']
	height_3 = isopter_labels_dict[skel_iso]['skeleton'][2][0]

	interpolated_height = (height_1*weight_1 + height_2*weight_2 + height_3*weight_3)

	return (interpolated_height)


				

def interpolate_three(isopter_labels_dict, min_distance_above, min_distance_below, third_dist, above_iso, below_iso, third_iso):
	inverse_dist_1 = 1/min_distance_above
	inverse_dist_2 = 1/min_distance_below
	inverse_dist_3 = 1/third_dist
	sum_inverse_dist = inverse_dist_1 + inverse_dist_2 + inverse_dist_3

	weight_1 = inverse_dist_1/sum_inverse_dist
	weight_2 = inverse_dist_2/sum_inverse_dist
	weight_3 = inverse_dist_3/sum_inverse_dist

	height_1 = isopter_labels_dict[above_iso]['height']
	height_2 = isopter_labels_dict[below_iso]['height']
	height_3 = isopter_labels_dict[third_iso]['height']

	interpolated_height = height_1*weight_1 + height_2*weight_2 + height_3*weight_3

	return (interpolated_height)


def find_distance_to_skeleton(isopter_labels_dict, circle_point, isopter):
	#print("DEBUG: find_distance_to_skeleton, CP = %s, isopter = %s" % (circle_point, isopter))
	min_distance = 10000
	min_point = []
	#print ("isopter: ",isopter_labels_dict[isopter])
	#print ("circle_point", circle_point)
	x_points = isopter_labels_dict[isopter]['skeleton'][0]
	y_points = isopter_labels_dict[isopter]['skeleton'][1]
	for i in range(0,len(x_points)):
		#print ("x_point, y_point")
		#print (x_points[i], y_points[i])
		current_distance = math.hypot(circle_point[0]-x_points[i],circle_point[1]-y_points[i])
		if (current_distance < min_distance):
			min_distance = current_distance
			min_point = [x_points[i], y_points[i]]
	if (min_distance <= 0):
		min_distance = 0.1
	#print ("min_distance, min_point")
	#print (min_distance, min_point)
	return (min_distance, min_point)

def find_distance_to_isopter(isopter_labels_dict, circle_point, isopter):
	#print("DEBUG: find_distance_to_isopter, CP = %s, isopter = %s" % (circle_point, isopter))
	min_distance = 10000
	min_point = []
	x_points = []
	y_points = []
	#print (isopter_labels_dict[isopter]['polygon'].exterior.coords)
	points = isopter_labels_dict[isopter]['polygon'].exterior.coords
	for point in points:
		x_points.append(point[0])
		y_points.append(point[1])
	for i in range(0,len(x_points)):
		current_distance = math.hypot(circle_point[0]-x_points[i],circle_point[1]-y_points[i])
		if (current_distance < min_distance):
			min_distance = current_distance
			min_point = [x_points[i], y_points[i]]
	if (min_distance <= 0):
		min_distance = 0.1
	return (min_distance, min_point)
	

def find_boundary(isopter_labels_dict, intersects_dict, isopter, other):
	below = ''
	in_child  = ''
	in_sc = ''
	isopter_children, scotoma_children = find_children(isopter_labels_dict, intersects_dict, isopter)
	for each in isopter_children:
		test_point = Point(other[0], other[1])
		polygon = isopter_labels_dict[each]['polygon']
		point_test = test_point.within(polygon)
		if (point_test == True):
			in_child = each
	for each in scotoma_children:
		test_point = Point(other[0], other[1])
		polygon = isopter_labels_dict[each]['polygon']
		point_test = test_point.within(polygon)
		if (point_test == True):
			in_sc = each
	if (in_child != ''):
		#print("DEBUG: inside(FB): ", isopter_labels_dict[in_child]['name'])
		return(find_boundary(isopter_labels_dict, intersects_dict, isopter_labels_dict[in_child], other))
	elif (in_sc != ''):
		#print("DEBUG: inside(FB): ", isopter_labels_dict[in_sc]['name'])
		return(find_boundary(isopter_labels_dict, intersects_dict, isopter_labels_dict[in_sc], other))
	else:
		min = 100000
		min_iso = ''
		#isopter_children, scotoma_children = find_children(isopter_labels_dict, intersects_dict, isopter)
		all_children = isopter_children + scotoma_children
		for each in all_children:
			x_points = []
			y_points = []
			points = isopter_labels_dict[each]['polygon'].exterior.coords
			for point in points:
				x_points.append(point[0])
				y_points.append(point[1])
			for i in range(0,len(x_points)):
				distance = math.hypot(other[0]-x_points[i],other[1]-y_points[i])
				if (distance < min):
					min = distance
					min_iso = each

		return (isopter, min_iso)

	
		
def find_children(isopter_labels_dict, intersects_dict, isopter):
	#if len(isopter) == 1:
	#	print("DEBUG: DEBUG: isopter = ", isopter_labels_dict[isopter]['name'])
	#else:
	#	print("DEBUG: DEBUG: isopter = ", isopter['name'])
	#print("DEBUG: Searching for children in ", isopter['name'])
	#print("DEBUG: Searching for children in ", isopter_labels_dict[isopter]['name'])
	iso_letter = ''
	## If we are calling the function with a letter/ID rather than the isopter "object" (I think)
	if (len(isopter) == 1):
		iso_letter = isopter
	else:
		if "letter" in isopter:
			iso_letter = isopter['letter']
		else:
			iso_letter = ''
			for each in isopter_labels_dict:
				isopter_labels_dict[each]['letter'] = each
				if (isopter['name'] == isopter_labels_dict[each]['name']):
					iso_letter = each
	isopter_children = intersects_dict[iso_letter]["isopters"]
	scotoma_children = intersects_dict[iso_letter]["scotomas"]
	#print("DEBUG:     FindChildren: I(%s), SC(%s)" % (isopter_children, scotoma_children))
	return (isopter_children, scotoma_children)

def combineIsopters(openIso, exIso):

#	print (openIso)
#	print (len(openIso))
#	print (exIso) 
	complete_isopters = []
	num_complete = 1
	open_len = len(openIso)
	extrap_len = len(exIso)

	if (open_len == extrap_len):
		for i in range(0,len(openIso)-3,4):
			open_name = openIso[i]
			#print ("open name")
			#print  (open_name)
			name = open_name.split("-")[0] + "-completed-" + str(num_complete)
			#print (name)
			#print ("open points")
			open_points = []
			for j in range(0,len(openIso[i+1])):
				open_points.append([openIso[i+1][j], openIso[i+2][j]])
			#	print ([openIso[i+1][j], openIso[i+2][j]])
			extrap_name = exIso[i]
			extrap_points = []
			#print ("extrapolated points")
			for j in range(0,len(exIso[i+1])):
				extrap_points.append([exIso[i+1][j], exIso[i+2][j]])
			#	print ([exIso[i+1][j], exIso[i+2][j]])
			open_object = LinearRing(open_points)
			extrap_object = LinearRing(extrap_points)
			if (open_object.is_ccw):
				open_object.coords = list(open_object.coords)[::-1]
			if (extrap_object.is_ccw):
				extrap_object.coords = list(extrap_object.coords)[::-1]
			#print (open_object.is_ccw)
			#print (extrap_object.is_ccw)
			#for point in open_object.coords[:]:
			#	print point
			complete_points = open_object.coords[:-1] + extrap_object.coords[:-1]
			#print ("complete_points")
			x_points = []
			y_points = []
			for point in complete_points:
				x_points.append(point[0])
				y_points.append(point[1])
			#	print (point[0], point[1])
			complete_isopters.append(name)
			complete_isopters.append(x_points)
			complete_isopters.append(y_points)
			complete_isopters.append([openIso[i+3][0]]*len(x_points))
	
	return (complete_isopters)
	
		


def parse_annotation(inputfile, outputfile, z_scale, debug):

	#print("in parseAnnotation")
	inputfile = inputfile
	#print(inputfile)
	dom = xml.dom.minidom.parse(inputfile)
	#print("test1")
	images = dom.getElementsByTagName("annotation-per-image")
	#print(images)
	for each in images:
		#Name of the output PDF
		fileName = get_imagename(each)
		if(debug == 1):
			print("FILE=", fileName)
		#Check for fiducials and find the scaling factor
		scale,num_fid,mX,mY,tX,tY,rX,rY,bX,bY,lX,lY,eye = fiducials_scale(each)
		if(debug == 1):
			print(scale,num_fid,mX,mY,tX,tY,rX,rY,bX,bY,lX,lY,eye)
		if (mX != 0 and mY != 0 and scale != 0):
			#Calculate average rotation and RRatio
			avgRotation, avgRRatio = rotation_calculation(mX,mY,tX,tY,rX,rY,bX,bY,lX,lY, debug)
			if(debug == 1):
				print("aveRotation=", avgRotation, " aveRRatio=", avgRRatio)
			#Rotate fiducial points
			rotatedFids, pts = rotate_fiducials(avgRotation,avgRRatio,mX,mY,tX,tY,rX,rY,bX,bY,lX,lY)


#def rotate_fiducials(avgRotation,avgRRatio,mX,mY,tX,tY,rX,rY,bX,bY,lX,lY):
#        fidX = [tX,rX,bX,lX]
#        fidY = [tY,rY,bY,lY]
#        fidX_pt = [90,0,270,180]
#        fidY_pt = [70,90,70,90]



			#Rotate isopter points and flip y-axis
			isopters, scotomas, openIso, exIso = collect_isopters(each,avgRotation,avgRRatio,mX,mY)
			print(isopters, flush=True)
			new_isos = combineIsopters(openIso, exIso)

			circle_points = []
			circle_degs = []
			#chunk_len = (float((mY-tY)/7))
			#degree_len = chunk_len/10
			degree_len = avgRRatio
			#Need 1-90 times degree length
			circle_distances = list(range(1,91))
			#print (circle_distances)
			circle_distances_deg = []
			for each in circle_distances:
				degree = each * degree_len
				circle_distances_deg.append(degree)

			circle_points.append([[mX,mY]])
			circle_degs.append([[0,0]])
			
			for each in circle_distances:
				circle_points.append(PointsInCircum(each,360,mX,mY))
				circle_degs.append(DegsInCircum(each,360))
	
			isopters = isopters + new_isos
			construct_hos(isopters, scotomas, circle_points, circle_degs, mX, mY, tX, tY, bX, bY, rX, rY, lX, lY, degree_len, outputfile, z_scale, eye)
			#print("Isopters = ", isopters)

def DegsInCircum(r,n):
	circle_degs = []

	for x in range(0,n+1):
		circle_degs.append([r,x])
	return (circle_degs)

def PointsInCircum(r,n,mX,mY):
	circle_points = []

	for x in range(0,n+1):
		circle_points.append([(math.cos(2*math.pi/n*x)*r)+mX,(math.sin(2*math.pi/n*x)*r)+mY])
		#print (((math.cos(2*math.pi/n*x)*r)+mX) , ((math.sin(2*math.pi/n*x)*r)+mY))
	return (circle_points)

def collect_isopters(each,avgRotation,avgRRatio,mX,mY):
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

	numI3eEx = 1
	numIII4eEx = 1
	numV4eEx = 1
	numIII4eExFull = 1
	numV4eExFull = 1

	for annotation in annotations:
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I1e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [120] * len(xlist2)
				zlist = [4] * len(xlist2)
				isopters.append('I1e ' + str(numI1e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				isopters.append(zlist)
				numI1e = numI1e + 1
			if (name == 'I1e-sc'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [120] * len(xlist2)
				zlist = [4] * len(xlist2)
				scotomas.append('I1e-sc' + str(numI1eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				scotomas.append(zlist)
				numI1eSC = numI1eSC + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I1e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [120] * len(xlist2)
				zlist = [4] * len(xlist2)
				openIso.append('I1e-open' + str(numI1eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				openIso.append(zlist)
				numI1eOpen = numI1eOpen + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I2e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [105] * len(xlist2)
				zlist = [3.5] * len(xlist2)
				isopters.append('I2e' + str(numI2e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				isopters.append(zlist)
				numI2e = numI2e + 1
			if (name == 'I2e-sc'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [105] * len(xlist2)
				zlist = [3.5] * len(xlist2)
				scotomas.append('I2e-sc' + str(numI2eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				scotomas.append(zlist)
				numI2eSC = numI2eSC + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I2e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [105] * len(xlist2)
				zlist = [3.5] * len(xlist2)
				openIso.append('I2e-open' + str(numI2eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				openIso.append(zlist)
				numI2eOpen = numI2eOpen + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I3e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [90] * len(xlist2)
				zlist = [3] * len(xlist2)
				isopters.append('I3e' + str(numI3e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				isopters.append(zlist)
				numI3e = numI3e + 1
			if (name == 'I3e-sc'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [90] * len(xlist2)
				zlist = [3] * len(xlist2)
				scotomas.append('I3e-sc' + str(numI3eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				scotomas.append(zlist)
				numI3eSC = numI3eSC + 1
			if (name == 'I3e-extrap'):
				xlist2, ylist2 = read_in_points_ex(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [90] * len(xlist2)
				zlist = [3] * len(xlist2)
				exIso.append('I3e-sc' + str(numI3eEx))
				exIso.append(xlist2)
				exIso.append(ylist2)
				exIso.append(zlist)
				numI3eEx = numI3eEx + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I3e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [90] * len(xlist2)
				zlist = [3] * len(xlist2)
				openIso.append('I3e-open' + str(numI3eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				openIso.append(zlist)
				numI3eOpen = numI3eOpen + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I4e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [75] * len(xlist2)
				zlist = [2.5] * len(xlist2)
				isopters.append('I4e' + str(numI4e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				isopters.append(zlist)
				numI4e = numI4e + 1
			if (name == 'I4e-sc'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [75] * len(xlist2)
				zlist = [2.5] * len(xlist2)
				scotomas.append('I4e-sc' + str(numI4eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				scotomas.append(zlist)
				numI4eSC = numI4eSC + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'I4e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [75] * len(xlist2)
				zlist = [2.5] * len(xlist2)
				openIso.append('I4e-open' + str(numI4eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				openIso.append(zlist)
				numI4eOpen = numI4eOpen + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'III4e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [45] * len(xlist2)
				zlist = [1.5] * len(xlist2)
				isopters.append('III4e' + str(numIII4e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				isopters.append(zlist)
				numIII4e = numIII4e + 1
			if (name == 'III4e-sc'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [45] * len(xlist2)
				zlist = [1.5] * len(xlist2)
				scotomas.append('III4e-sc' + str(numIII4eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				scotomas.append(zlist)
				numIII4eSC = numIII4eSC + 1
			if (name == 'III4e-extrap'):
				xlist2, ylist2 = read_in_points_ex(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [45] * len(xlist2)
				zlist = [1.5] * len(xlist2)
				exIso.append('III4e-extrap' + str(numIII4eEx))
				exIso.append(xlist2)
				exIso.append(ylist2)
				exIso.append(zlist)
				numIII4eEx = numIII4eEx + 1
			if (name == 'III4e-extrap-full'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [45] * len(xlist2)
				zlist = [1.5] * len(xlist2)
				isopters.append('III4e-extrap-full' + str(numIII4eExFull))
				isopters.append(xlist2)
				isopters.append(ylist2)
				isopters.append(zlist)
				numIII4eExFull = numIII4eExFull + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'III4e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [45] * len(xlist2)
				zlist = [1.5] * len(xlist2)
				openIso.append('III4e-open' + str(numIII4eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				openIso.append(zlist)
				numIII4eOpen = numIII4eOpen + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "SplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'V4e'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [15] * len(xlist2)
				zlist = [0.5] * len(xlist2)
				isopters.append('V4e' + str(numV4e))
				isopters.append(xlist2)
				isopters.append(ylist2)
				isopters.append(zlist)
				numV4e = numV4e + 1
			if (name == 'V4e-sc'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [15] * len(xlist2)
				zlist = [0.5] * len(xlist2)
				scotomas.append('V4e-sc' + str(numV4eSC))
				scotomas.append(xlist2)
				scotomas.append(ylist2)
				scotomas.append(zlist)
				numV4eSC = numV4eSC + 1
			if (name == 'V4e-extrap'):
				xlist2, ylist2 = read_in_points_ex(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [15] * len(xlist2)
				zlist = [0.5] * len(xlist2)
				exIso.append('V4e-extrap' + str(numV4eEx))
				exIso.append(xlist2)
				exIso.append(ylist2)
				exIso.append(zlist)
				numV4eEx = numV4eEx + 1
			if (name == 'V4e-extrap-full'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [15] * len(xlist2)
				zlist = [0.5] * len(xlist2)
				isopters.append('V4e-extrap-full' + str(numV4eExFull))
				isopters.append(xlist2)
				isopters.append(ylist2)
				isopters.append(zlist)
				numV4eExFull = numV4eExFull + 1
		if (annotation.getElementsByTagName("controltype")[0].firstChild.nodeValue == "PolysplineRoi"):
			name = annotation.getElementsByTagName("name")[0].firstChild.nodeValue
			name = name.replace(" ","")
			if (name == 'V4e-open'):
				xlist2, ylist2 = read_in_points(annotation,avgRotation,avgRRatio,mX, mY)
				#zlist = [15] * len(xlist2)
				zlist = [0.5] * len(xlist2)
				openIso.append('V4e-open' + str(numV4eOpen))
				openIso.append(xlist2)
				openIso.append(ylist2)
				openIso.append(zlist)
				numV4eOpen = numV4eOpen + 1

	return (isopters, scotomas, openIso, exIso)

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
	#print("XLIST: ", xlist2)
	#print("YLIST: ", ylist2)
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
	for i in range(len(ylist2)):
		ylist2[i] = (ylist2[i] - mY)*(-1)
		ylist2[i] = ylist2[i] + mY
	return (xlist2, ylist2)

def rotate_points(avgRotation,avgRRatio,mX, mY, xlist, ylist):
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

def rotate_fiducials(avgRotation,avgRRatio,mX,mY,tX,tY,rX,rY,bX,bY,lX,lY):
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
	#for i in range(len(rFidX)):
		fidComplete.append(xTest)
		fidComplete.append(yTest)
		#fidComplete.append(rFidX[i])
		#fidComplete.append(rFidY[i])
		ptsComplete.append(rFidX_pt[i])
		ptsComplete.append(rFidY_pt[i])
	print("DEBUG: tX=",tX," tY=",tY, flush=True)
	print("DEBUG: fidC=",fidComplete," ptsComplete=",ptsComplete, flush=True)
	return (fidComplete, ptsComplete)

def rotation_calculation(mX,mY,tX,tY,rX,rY,bX,bY,lX,lY, debug):
	length = 0
	avgRotation = 0
	avgRRatio = 0
	dist_fid = [tX,tY,rX,rY,bX,bY,lX,lY]
	dist_loc = [90,70,0,90,270,70,180,90]
	fid = []
	loc = []
	angles = []

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
		cx = fid[i] - mX
		cy = mY - fid[i+1]
		#pixel polar coordinates
		ptheta = math.atan2((cy),cx)
		pr = math.hypot(cx,cy)
		#GVF coordinate conversion set up
		thetaDifference = (ptheta - (math.radians(loc[i])))
		angles.append(thetaDifference)
		#while (abs(thetaDifference)>math.pi):
		#	if(thetaDifference<0):
		#		thetaDifference += 2*math.pi
		#	else:
		#		thetaDifference -= 2*math.pi
		if(debug == 1):
			print("thetaDifference=", thetaDifference)
		#avgRotation += (thetaDifference/length)
		avgRRatio += (pr/loc[i+1])/length
	avgRotation = average_angles(angles)
	if(debug == 1):
		print("avgRotation=", avgRotation)
	while (abs(avgRotation)>math.pi):
		if(avgRotation<0):
			avgRotation += 2*math.pi
		else:
			avgRotation -= 2*math.pi
	return (avgRotation, avgRRatio)

def average_angles(angles):
	sum_cos = 0;
	sum_sin = 0
	for i in range(0,len(angles)):
		sum_cos += math.cos(angles[i])
		sum_sin += math.sin(angles[i])
	return(math.atan2(sum_sin, sum_cos))

	
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
	eye = ''

	annotations = each.getElementsByTagName("annotation")
	for annotation in annotations:
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue.startswith("90,70")):
			topX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			topY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue.startswith("270,70")):
			bottomX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			bottomY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue.startswith("180,90")):
			leftX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			leftY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue.startswith("0,90")):
			rightX = int(annotation.getElementsByTagName("x")[0].firstChild.nodeValue)
			rightY = int(annotation.getElementsByTagName("y")[0].firstChild.nodeValue)
			num_fid = num_fid + 1
		if(annotation.getElementsByTagName("name")[0].firstChild.nodeValue.startswith("0,0")):
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
			distance = (middleY - topY)
		elif (bottomX != 0 and bottomY != 0):
			distance = (bottomY - middleY)
		elif (leftX != 0 and leftY != 0):
			distance = ((7/9)*(middleX - leftX))
		elif (rightX != 0 and rightY != 0):
			distance = ((7/9)*(rightX - middleX))
		else:
			distance = 0
	## Needed for images that are upside-down
	distance = abs(distance)
	if (distance > 0):
		size = distance/300
		scale = 25*size
	
	return(scale,num_fid,middleX,middleY,topX,topY,rightX,rightY,bottomX,bottomY,leftX,leftY,eye)

def get_imagename(each):
	imagename = each.getElementsByTagName("imagename")
	oldFileName = imagename[0].firstChild.nodeValue
	newFileName = []
	oldFileName = oldFileName.split('/')
	if len(oldFileName) > 1:
		oldFileName = oldFileName[1]
	else:
		oldFileName = oldFileName[0]
	newFileName = ''.join(oldFileName)
	newFileName = newFileName.split('_')
	if len(newFileName) > 1:
		newFileName = newFileName[1]
	else:
		newFileName = newFileName[0]
	newFileName = ''.join(newFileName)
	fileName = oldFileName
	return(fileName)

def main(argv):
	inputfile = ''
	outputfile = ''
	debug = 0
	## May need to set a better default than 1.0 for perception
	z_scale = 1.0
	try:
		#opts, args = getopt.getopt(argv,"huo:pi:s",["input=", "uiowa-default-colors","output-base-name=", "pdf-only"])
		opts, args = getopt.getopt(argv,"hdo:i:z:",["input=", "output=", "scale=" ])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			usage()
			sys.exit()
		elif opt == '-d':
			debug = 1
		elif opt in ("-i", "--input"):
			inputfile = arg
			## set default outputfile
			tmpstring = inputfile + ".html"
			outputfile = re.sub(".txt.html$", ".html", tmpstring)
			#outputfile = inputfile + ".html"
			print ('Input file is ', inputfile)
		elif opt in ("-o", "--output"):
			outputfile = arg
		elif opt in ("-z", "--scale"):
			z_scale = arg
	if(inputfile == ''):
		usage()
		sys.exit()
	try:
		z_scale = float(z_scale)
	except:
		print("Error: Z-scale must be numeric! (",z_scale," specified)")
		usage()
		sys.exit()
	print ('Output file is ', outputfile)
	#print ("Z scale = ", z_scale)

	parse_annotation(inputfile, outputfile, z_scale, debug)

def usage():
	print (__doc__)
if __name__ == "__main__":
	main(sys.argv[1:])
