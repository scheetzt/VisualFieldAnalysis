
import argparse
import os
import glob
import re
from pprint import pprint
import xml.etree.ElementTree as ET
from tkinter import *
from PIL import Image, ImageDraw, ImageColor, ImageTk
#import ImageTk
from math import atan2, degrees

## set to defaults
rotation = 0
scale = 1.0
offset_x = 0
offset_y = -70
QA_pass = "NA"
fiducial_pass = "NA"
VF_pass = "NA"
partial_field = "NA"
##declare global data
image_dict = dict()
render_dict = dict()
xml_dict = dict()
all_dict = dict()
current_image = 0
show_image = 0
final_list = list()
num_images = 0
per_image_data = dict()
fast = 1
## create the Tk Window
TkRoot = Tk()
canvas = Canvas(TkRoot, width=1200, height=900, bg='white')
canvas.pack(expand=YES, fill=BOTH)
outputfile = open("output.txt", "w")

class MyFP:
	def __init__(self, type, x, y):
		self.type = type
		self.x = x
		self.y = y
	def getX(self):
		return(self.x)
	def getY(self):
		return(self.y)

def onKeyPress(event):
	global rotation, scale, show_image, current_image, final_list
	global im2, im1, image, canvas, num_images, outputfile
	global offset_x, offset_y, per_image_data
	global QA_pass, fiducial_pass, VF_pass, partial_field
	global fast
	if(event.char == "?"):
		print("Usage:")
		print("  QA: [g] = good, [b] = bad")
		print("  VF quality: [a] = good, [z] = bad")
		print("  VF Type: [i] = partial/incomplete")
		print("  Fiducial quality: [f] = good")
		print("  Controls: [q] = quit, [w] = write")
		print("  Movement: [h|e] = left, [l|r] = right, [j] = down, [k] = up")
		print("  Rotation: [c] = clockwise, [v] = counterclockwise")
		print("  Scale: [d] = toggles fast mode on/off")
		print("  Scale: [u] = bigger, [y] = smaller")
		print("  Stack: [n] = next, [p] = previous")
		print("  Overlay: [s] = toggles showing the rendered VF")
	if(event.char == "g"):
		outputfile.write("QA: good " + final_list[current_image])
		QA_pass = "Y"
		print("QA: good " + final_list[current_image]+"  QA="+QA_pass)
	if(event.char == "b"):
		outputfile.write("QA: bad " + final_list[current_image])
		QA_pass = "N"
		print("QA: bad " + final_list[current_image])
	if(event.char == "z"):
		outputfile.write("QA: VF bad " + final_list[current_image])
		VF_pass = "N"
		print("QA: VF bad " + final_list[current_image])
	if(event.char == "a"):
		outputfile.write("QA: VF good " + final_list[current_image])
		VF_pass = "Y"
		print("QA: VF good " + final_list[current_image])
	if(event.char == "i"):
		if(partial_field == "Y"):
			partial_field = "N"
			outputfile.write("QA: Not a Partial Field " + final_list[current_image])
			print("QA: Not a Partial Field " + final_list[current_image])
		else:
			partial_field = "Y"
			outputfile.write("QA: Partial Field " + final_list[current_image])
			print("QA: Partial Field " + final_list[current_image])
	if(event.char == "f"):
		if (fiducial_pass == "Y"):
			outputfile.write("QA: fiducials bad " + final_list[current_image])
			fiducial_pass = "N"
			print("QA: fiducials bad " + final_list[current_image])
		else:	
			outputfile.write("QA: fiducials good " + final_list[current_image])
			fiducial_pass = "Y"
			print("QA: fiducials good " + final_list[current_image])
	if(event.char == "q"):
		print("entered: q. Exiting...")
		writeQAPerImage()
		exit()
	if(event.char == "w"):
		print("entered: w. Writing per_image_data...")
		writeQAPerImage()
	if(event.char == "d"):
		if(fast == 1):
			fast = 10
		else:
			fast = 1
		print("entered: d (fast)")
		print(fast)
	if(event.char == "c"):
		rotation += 1
		print("entered: c")
		print(rotation)
	if(event.char == "v"):
		rotation -= 1
		print("entered: v")
		print(rotation)
	if(event.char == "u"):
		scale += 0.01 * fast
		print("entered: u")
		print(scale)
	if(event.char == "y"):
		scale -= 0.01 * fast
		print("entered: o")
		print(scale)
	if(event.char == "r" or event.char == "l"):
		offset_x += 1 * fast
		print("entered: right")
		print(offset_x)
	if(event.char == "e" or event.char == "h"):
		offset_x -= 1 * fast
		print("entered: left")
		print(offset_x)
	if(event.char == "j"):
		offset_y += 1 * fast
		print("entered: j")
		print(offset_y)
	if(event.char == "k"):
		offset_y -= 1 * fast
		print("entered: k")
		print(offset_y)
	if(event.char == "n"):
		print("entered: n")
		setPerImage()
		## increment the counter if not at end
		if(current_image < num_images - 1):
			current_image += 1
		else:
			print("at end already")
		## reset to default rotation, retain previous scaling (often correct)
		rotation = 0
		QA_pass = "NA"
		VF_pass = "NA"
		fiducial_pass = "NA"
		partial_field = "NA"
		## now load values for the next image
		getPerImage()
		print("" + str(current_image) + ": " + final_list[current_image])
	if(event.char == "p"):
		print("entered: p")
		setPerImage()
		if(current_image > 0):
			current_image -= 1
		else:
			print("at beginning already")
		## reset to default rotation, retain previous scaling (often correct)
		rotation = 0
		QA_pass = "NA"
		VF_pass = "NA"
		fiducial_pass = "NA"
		partial_field = "NA"
		## now load values for the next image
		getPerImage()
		print("" + str(current_image) + ": " + final_list[current_image])
	if(event.char == "s"):
		if(show_image == 1):
			show_image = 0
		else:
			show_image = 1
		print("entered: s")
		print(show_image)

	renderImage()

########################################################################
def setPerImage():
	global final_list, current_image, per_image_data
	global offset_x, offset_y, scale, rotation
	global fiducial_pass, QA_pass, VF_pass, partial_field

	imagename = final_list[current_image]
	if imagename not in per_image_data:
		per_image_data[imagename] = dict()
	## save the relevant data
	per_image_data[imagename].update({"rotation" : rotation})
	per_image_data[imagename].update({"scale" : scale})
	per_image_data[imagename].update({"offset_x" : offset_x})
	per_image_data[imagename].update({"offset_y" : offset_y})
	per_image_data[imagename].update({"fiducial_pass" : fiducial_pass})
	per_image_data[imagename].update({"QA_pass" : QA_pass})
	per_image_data[imagename].update({"VF_pass" : VF_pass})
	per_image_data[imagename].update({"partial_field" : partial_field})
	#print("IM: "+imagename+"  QA: "+QA_pass)
	#pprint(per_image_data[imagename])

########################################################################
def getPerImage():
	global final_list, current_image, per_image_data
	global offset_x, offset_y, scale, rotation, field_scale
	global fiducial_pass, QA_pass, xml_dict, VF_pass, partial_field

	imagename = final_list[current_image]
	if imagename in per_image_data:
		if "rotation" in per_image_data[imagename]:
			rotation = per_image_data[imagename]["rotation"]
		if "scale" in per_image_data[imagename]:
			scale = per_image_data[imagename]["scale"]
		if "offset_x" in per_image_data[imagename]:
			offset_x = per_image_data[imagename]["offset_x"]
		if "offset_y" in per_image_data[imagename]:
			offset_y = per_image_data[imagename]["offset_y"]
		if "fiducial_pass" in per_image_data[imagename]:
			fiducial_pass = per_image_data[imagename]["fiducial_pass"]
		if "QA_pass" in per_image_data[imagename]:
			QA_pass = per_image_data[imagename]["QA_pass"]
		if "VF_pass" in per_image_data[imagename]:
			VF_pass = per_image_data[imagename]["VF_pass"]
		if "partial_field" in per_image_data[imagename]:
			partial_field = per_image_data[imagename]["partial_field"]
	else:
		#print("Inside getPerImage() else.")
		## Need the (0,0), (90,70), and (270,70) fiducial points
		if "0,0" in xml_dict[imagename] and "90,70" in xml_dict[imagename] and "270,70" in xml_dict[imagename]:
			#print("Inside getPerImage() else, if.")
			## first determine scale
			y1 = xml_dict[imagename]["90,70"].getY()
			x1 = xml_dict[imagename]["90,70"].getX()
			y2 = xml_dict[imagename]["270,70"].getY()
			x2 = xml_dict[imagename]["270,70"].getX()
			height = abs(y1-y2)
			scale = height / 681.0

			## compute field_scale
			temp_image = Image.open(image_dict[imagename])
			field_scale = 1.0
			if(temp_image.width > 1200):
				field_scale = 1200/temp_image.width
			if(temp_image.height > 900):
				if(900/temp_image.height < field_scale):
					field_scale = 900/temp_image.height

			## then determine offset_x and offset_y
			y0 = xml_dict[imagename]["0,0"].getY()*field_scale
			x0 = xml_dict[imagename]["0,0"].getX()*field_scale

			offset_x = int(x0 - (487*scale*field_scale))
			offset_y = int(y0 - (420*scale*field_scale))
			## then determine rotation
			#rotation = math.degrees(math.atan2(x2-x1, y2-y1))
			rotation = degrees(atan2(x2-x1, y2-y1))


########################################################################
## read data in from QA file
## FILENAME  fiducial_flag  good_flag  offset_x  offset_y  scale  rotation
def readQAPerImage():
	global per_image_data

	qafile = open("file.QAdata_v15.txt", "r")
	for line in qafile:
		line.rstrip()		# remove trailing white-space
		arr = re.split(r'\t', line)
		imagename = arr[0]
		if len(arr) != 9:       ## updated for v15 format, which includes 9 fields
			print("ERROR: Bad line: should have 9 fields, found "+str(len(arr))+"\n"+line)
		elif imagename == "Imagename":
			print("WARNING: Skipping header line.")
		else:
			if imagename not in per_image_data:
				per_image_data[imagename] = dict()
			per_image_data[imagename].update({"rotation" : int(arr[8])})
			per_image_data[imagename].update({"scale" : float(arr[7])})
			per_image_data[imagename].update({"offset_y" : int(arr[6])})
			per_image_data[imagename].update({"offset_x" : int(arr[5])})
			per_image_data[imagename].update({"partial_field" : arr[4]})
			per_image_data[imagename].update({"VF_pass" : arr[3]})
			per_image_data[imagename].update({"QA_pass" : arr[2]})
			per_image_data[imagename].update({"fiducial_pass" : arr[1]})

########################################################################
def writeQAPerImage():
	global per_image_data
	qafile = open("file.QAdata_v15.txt", "w")
	qafile.write("Imagename\tFiducial_pass\tQA_pass\tVF_pass\tpartial\toffset_x\toffset_y\tscale\trotation\n")
	for imagename in per_image_data:
		#print("IM: "+imagename)
		#pprint(per_image_data[imagename])
		string = "" + imagename 
		if "fiducial_pass" in per_image_data[imagename]:
			string = string + "\t" +  per_image_data[imagename]["fiducial_pass"]
		else:
			string = string + "\tNA"
		if "QA_pass" in per_image_data[imagename]:
			string = string + "\t" +  per_image_data[imagename]["QA_pass"]
		else:
			string = string + "\tNA"
		if "VF_pass" in per_image_data[imagename]:
			string = string + "\t" +  per_image_data[imagename]["VF_pass"]
		else:
			string = string + "\tNA"
		if "partial_field" in per_image_data[imagename]:
			string = string + "\t" +  per_image_data[imagename]["partial_field"]
		else:
			string = string + "\tNA"
		if "offset_x" in per_image_data[imagename]:
			string = string + "\t" +  str(per_image_data[imagename]["offset_x"])
		else:
			string = string + "\tNA"
		if "offset_y" in per_image_data[imagename]:
			string = string + "\t" +  str(per_image_data[imagename]["offset_y"])
		else:
			string = string + "\tNA"
		if "scale" in per_image_data[imagename]:
			string = string + "\t" +  str(per_image_data[imagename]["scale"])
		else:
			string = string + "\tNA"
		if "rotation" in per_image_data[imagename]:
			string = string + "\t" +  str(int(per_image_data[imagename]["rotation"]))
		else:
			string = string + "\tNA"
		string = string + "\n"
		qafile.write(string)
	qafile.close()



########################################################################
def parseArguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('--imageDir')
	parser.add_argument('--renderDir')
	parser.add_argument('--xmlFile')
	args = parser.parse_args()

	print(args)
	arg_flag = 0
	if(args.imageDir is None):
		#print("Error: You must specify imageDir")
		arg_flag = 1
	if(args.renderDir is None):
		#print("Error: You must specify renderDir")
		arg_flag = 1
	if(args.xmlFile is None):
		#print("Error: You must specify xmlFile")
		arg_flag = 1
	if(arg_flag == 1):
		parser.print_help()
		exit()
	return(args.imageDir, args.renderDir, args.xmlFile)

def checkArguments(image_path, render_path, xmlFile):
	## verify that imageDir is a directory, and contains images
	if(os.path.exists(image_path) is not True or os.path.isdir(image_path) is not True):
		print("Error: the specified imageDir is not a directory, or does not exist")
	if(os.path.exists(render_path) is not True or os.path.isdir(render_path) is not True):
		print("Error: the specified renderDir is not a directory, or does not exist")
	if(os.path.exists(xmlFile) is not True):
		print("Error: the specified xmlFile does not exist")

def checkFiles(image_list, render_list, xml_file):
	image_dict = dict()
	all_dict = dict()
	for x in image_list:
		y = os.path.basename(x)
		image_dict[y] = x
		y = re.sub(".JPG", "", y)
		y = re.sub(".PNG", "", y)
		all_dict[y] = 1
	render_dict = dict()
	for x in render_list:
		y = os.path.basename(x)
		## convert to PNG
		png_file = x
		png_file = re.sub(".pdf", ".png", png_file)
		if(os.path.exists(png_file) is not True):
			CMD = "convert " + x + " -background white -flatten -alpha set "
			CMD = CMD + " -channel A -evaluate set 50% " + png_file
			os.system(CMD)
		render_dict[y] = png_file
		y = re.sub(".JPG.pdf", "", y)
		all_dict[y] = 1
	pprint(image_dict)
	pprint(render_dict)

	xml_dict = dict()
	tree = ET.parse(xml_file)
	root=tree.getroot()
	for child in root:
		imagename = ""
		if child.tag == "annotation-per-image":
			#print("MAIN", child.tag, child.text)
			# find the image name, and the fiducial points
			for child2 in child:
				if child2.tag == "imagename":
					imagename = child2.text
					imagename = os.path.basename(imagename)
					#imagename = re.sub(".JPG", "", imagename)
					#print("NAME = ", imagename)
					## create a new dict entry
					xml_dict[imagename] = dict()
					#f_points[imagename] = dict()
				if child2.tag == "annotations":
					for annot_child in child2:
						node = annot_child.find("name")
						node_name = node.text
						node_name = node_name.strip()
						#print("NODE = (", node_name, ")")
						if node_name == "0,0" or node_name == "0,90" or node_name == "90,70" or node_name == "180,90" or node_name == "270,70":
							roi_node = annot_child.find("roi")
							mark_node = roi_node.find("mark")
							x_node = mark_node.find("x")
							y_node = mark_node.find("y")
							x = int(x_node.text)
							y = int(y_node.text)
							#print("\t%s\tX=%d\tY=%d" % (node_name, x, y))
							#print("\t", node_name, "X=", x, "Y=", y)
							tmp = MyFP(node_name, x, y)
							xml_dict[imagename].update({node_name : tmp})
	
	match_count = 0
	keys_to_remove = dict()
	print("DEBUG: before key removal")
	## build dict of all images in source_dir and/or render_dir
	for x in all_dict:
		yj = x + ".JPG"
		yp = x + ".PNG"
		z = x + ".pdf"
		y = ""
		if yj in image_dict:
			y = yj
			z = x + ".JPG.pdf"
		elif yp in image_dict:
			y = yp
			z = x + ".PNG.pdf"
		else:
			keys_to_remove.update({x:1})
			print("Removing key. Not in image_dict. " + x)
		## remove unless both image and render are available
		if y not in image_dict or z not in render_dict:
			keys_to_remove.update({x:1})
			print("Removing key. Not in render_dict. " + x)
		## also remove if no Fiducial points
		if y in image_dict:
			#z = image_dict[y]
			print("DEBUG: checking xml for " + y)
			if y in xml_dict:
				ttt = len(xml_dict[y])
				print("# FP = " + str(ttt))
				if ttt < 1:
					keys_to_remove.update({x:1})
					print("Removing key. Not enough FPs. " + y)
			else:
				print(".. did not find "+z+" in xml_dict")
				keys_to_remove.update({x:1})
				print("Removing key. No XML. " + y)
	for x in keys_to_remove:
		print("..removing key: "+x)
		del all_dict[x]
	if(len(all_dict) < 1):
		print("ERROR: No matching source and rendered images available\nExiting...\n")
		exit()
	return(image_dict, render_dict, xml_dict, all_dict)

def renderImage():
	#global canvas
	global image_dict, render_dict, xml_dict, final_list
	global rotation, scale, show_image
	global field_scale
	global im2, im1, image, canvas
	global offset_x, offset_y, per_image_data
	#pprint(final_list[current_image])

	## Create the PIL image (im1) for initial drawing/rendering
	image_file = image_dict[final_list[current_image]]
	#print("Image file 1 = " + image_file)	
	im1 = Image.open(image_file)

	## scale im1 if larger than 1200 x 900 (which is the size of the view)
	field_scale = 1.0
	if(im1.width > 1200):
		field_scale = 1200/im1.width
	if(im1.height > 900):
		if(900/im1.height < field_scale):
			field_scale = 900/im1.height
	if(field_scale < 1.0):
		height = int(im1.height * field_scale)
		width = int(im1.width * field_scale)
		im1 = im1.resize(size=(width,height))

	temp = final_list[current_image]
	temp = re.sub(".JPG", ".JPG.pdf", temp)
	temp = re.sub(".PNG", ".PNG.pdf", temp)

	#print("Image file 2 = " + render_dict[temp])	
	im2 = Image.open(render_dict[temp])
	height = int(im2.height * scale * field_scale)
	width = int(im2.width * scale * field_scale)
	im2 = im2.resize(size=(width,height))
	im2 = im2.rotate(rotation)
	if(show_image == 1):
		im1.paste(im2, (offset_x,offset_y), im2.convert("RGBA"))

	## get Fiducial points
	#pprint(xml_dict)
	#if image_file in xml_dict:
	#print(final_list[current_image])
	if final_list[current_image] in xml_dict:
		#print("Fiducials found!\n")
		draw = ImageDraw.Draw(im1)
		for fp in xml_dict[final_list[current_image]]:
			name = fp
			x = xml_dict[final_list[current_image]][fp].getX()
			y = xml_dict[final_list[current_image]][fp].getY()
			#print("NODE: " + fp +" "+ str(x) +" "+ str(y))
			if(field_scale < 1.0):
				x = x * field_scale
				y = y * field_scale
			drawFiducial(im1, name, x, y)
	else:
		print("No fiducials.\n")

	## Draw the 0,0
	#draw.line((526, 393, 546, 413), width=3, fill=(255,0,0))
	#draw.line((526, 413, 546, 393), width=3, fill=(255,0,0))

	## create the ImageTk compliant version of im1, bind it to the canvas
	image = ImageTk.PhotoImage(im1)
	canvas.create_image(0, 0, image=image, anchor=NW)

def drawFiducial(image, name, real_x, real_y):
	draw = ImageDraw.Draw(image)
	start_x = real_x
	start_y = real_y
	#print("FP(%s, %s) = %d, %d" % (image, name, start_x, start_y))
	end_x = start_x
	end_y = start_y
	if name == "0,0":
		start_x = real_x - 10
		end_x = real_x + 10
		start_y = real_y - 10
		end_y = real_y + 10
		draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
		start_x = real_x - 10
		end_x = real_x + 10
		start_y = real_y + 10
		end_y = real_y - 10
		draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
	if name == "0,90":
		start_x = real_x
		end_x = real_x + 10
		start_y = real_y
		end_y = real_y + 10
		draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
		start_x = real_x
		end_x = real_x + 10
		start_y = real_y
		end_y = real_y - 10
		draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
	if name == "90,70":
		start_x = real_x
		end_x = real_x + 10
		start_y = real_y
		end_y = real_y - 10
		draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
		start_x = real_x
		end_x = real_x - 10
		start_y = real_y
		end_y = real_y - 10
		draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
	if name == "180,90":
		start_x = real_x
		end_x = real_x - 10
		start_y = real_y
		end_y = real_y + 10
		draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
		start_x = real_x
		end_x = real_x - 10
		start_y = real_y
		end_y = real_y - 10
		draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
	if name == "270,70":
		start_x = real_x
		end_x = real_x - 10
		start_y = real_y
		end_y = real_y + 10
		draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
		start_x = real_x
		end_x = real_x + 10
		start_y = real_y
		end_y = real_y + 10
		draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
	del draw


def main():
	global image_dict, render_dict, xml_dict, final_list, canvas
	global num_images
	image_path = ""
	render_path = ""
	xmlFile = ""
	(image_path, render_path, xmlFile) = parseArguments()

	checkArguments(image_path, render_path, xmlFile)

	outputfile.write("test")

	## build a list of the files in the image_path and render_path
	temp_string = image_path + '/*.JPG'
	image_list_jpg = glob.glob(temp_string)
	temp_string = image_path + '/*.PNG'
	image_list_png = glob.glob(temp_string)
	print("# JPG = " +str(len(image_list_jpg)))
	print("# PNG = " +str(len(image_list_png)))
	image_list = image_list_jpg + image_list_png
	if(len(image_list) == 0):
		temp_string = image_path + '/*/*.JPG'
		image_list_jpg = glob.glob(temp_string)
		temp_string = image_path + '/*/*.PNG'
		image_list_png = glob.glob(temp_string)
		image_list = image_list_jpg + image_list_png
		print("# JPG = " +str(len(image_list_jpg)))
		print("# PNG = " +str(len(image_list_png)))
		if(len(image_list)==0):
			print("ERROR: no images found!\n")
			exit()
	#pprint(image_list)
	temp_string = render_path + '/*.pdf'
	render_list = glob.glob(temp_string)
	if(len(render_list) == 0):
		temp_string = render_path + '/*/*.pdf'
		render_list = glob.glob(temp_string)
		if(len(render_list)==0):
			print("ERROR: no renders found!\n")
			exit()
	#pprint(render_list)

	## verify we have matched pairs of images and rendered traces
	(image_dict, render_dict, xml_dict, all_dict) = checkFiles(image_list, render_list, xmlFile)
	print("Success. After checkFiles().\n")

	for x in sorted(all_dict):
		xj = x + ".JPG"
		xp = x + ".PNG"
		if xj in image_dict:
			final_list.append(xj)
		elif xp in image_dict:
			final_list.append(xp)


	## read in QA per-image-file
	if os.path.isfile("file.QAdata_v15.txt"):
		readQAPerImage()

	num_images = len(final_list)
	# load parameters, if available, for the first image
	getPerImage()

	image_file = image_dict[final_list[current_image]]
	#print("Image file = " + image_file)	
	im1 = Image.open(image_file)
	image = ImageTk.PhotoImage(im1)
	#pprint(canvas)
	canvas.create_image(0, 0, image=image, anchor=NW)
	renderImage()

	## bind onKeyPress() to respond to keyboard events
	TkRoot.bind('<KeyPress>', onKeyPress)
	TkRoot.focus()
	mainloop()




if __name__ == '__main__':
	main()










#for image in f_points.keys():
#	im = Image.open(source_dict[os.path.basename(image)])
#	draw = ImageDraw.Draw(im)
#	for fp in f_points[image].keys():
#		#print("FP(%s, %s)" % (image, fp))
#		real_x = f_points[image][fp].getX()
#		real_y = f_points[image][fp].getY()
#		start_x = real_x
#		start_y = real_y
#		print("FP(%s, %s) = %d, %d" % (image, fp, start_x, start_y))
#		end_x = start_x
#		end_y = start_y
#		if fp == "0,0":
#			start_x = real_x - 10
#			end_x = real_x + 10
#			start_y = real_y - 10
#			end_y = real_y + 10
#			draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
#			start_x = real_x - 10
#			end_x = real_x + 10
#			start_y = real_y + 10
#			end_y = real_y - 10
#			draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
#		if fp == "0,90":
#			start_x = real_x 
#			end_x = real_x + 10
#			start_y = real_y 
#			end_y = real_y + 10
#			draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
#			start_x = real_x 
#			end_x = real_x + 10
#			start_y = real_y 
#			end_y = real_y - 10
#			draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
#		if fp == "90,70":
#			start_x = real_x 
#			end_x = real_x + 10
#			start_y = real_y 
#			end_y = real_y - 10
#			draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
#			start_x = real_x 
#			end_x = real_x - 10
#			start_y = real_y 
#			end_y = real_y - 10
#			draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
#		if fp == "180,90":
#			start_x = real_x 
#			end_x = real_x - 10
#			start_y = real_y 
#			end_y = real_y + 10
#			draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
#			start_x = real_x 
#			end_x = real_x - 10
#			start_y = real_y 
#			end_y = real_y - 10
#			draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
#		if fp == "270,70":
#			start_x = real_x 
#			end_x = real_x - 10
#			start_y = real_y 
#			end_y = real_y + 10
#			draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
#			start_x = real_x 
#			end_x = real_x + 10
#			start_y = real_y 
#			end_y = real_y + 10
#			draw.line((start_x, start_y, end_x, end_y), width=3, fill=(255,0,0))
#
#
#
#
#
#
#
