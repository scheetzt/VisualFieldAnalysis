#!/usr/bin/python

## Converts Truthmarker generated XML file into individual files, allowing them to be fixed
## Usage: python3 extractAnnotationFiles.py [options] [source]
## Options:
##   -i ..., --input=...             Use specified XML file.
##   -h                              print out a help message describing the options
##   -o DIR, --output-dir DIR        specify the directory in which output will be rendered

## DOC STRING attempt
"""Usage: extractAnnotationFiles.py [-h] [-i FILE | --input FILE] [-o DIR | --output-dir DIR] """


import os, sys, getopt, pprint, re, math
import xml.dom.minidom
import os
import datetime
import string
from pprint import pprint



#######


def get_imagepath(each, output_dir):
	imagename = each.getElementsByTagName("imagename")
	fileName = imagename[0].firstChild.nodeValue
	fileName = fileName.split('/')
	if len(fileName) > 1:
		fileName = fileName[-1]
	else:	
		fileName = fileName[0]
	return(output_dir + "/" + fileName + ".xml")

## the real meat...

def extractAnnotationFiles(inputfile, output_dir):
	inputfile = inputfile
	dom = xml.dom.minidom.parse(inputfile)
	images = dom.getElementsByTagName("annotation-per-image")

	for each in images:
		#Name of the output PDF
		fileName = get_imagepath(each, output_dir)
		output = open(fileName, "w")
		pretty_xml_as_string = each.toprettyxml()
		pretty_xml_as_string =  os.linesep.join([s for s in pretty_xml_as_string.splitlines() if s.strip()])
		# Add <?xml version="1.0"?>, and a newline at the end
		pretty_xml_as_string = "<?xml version=\"1.0\"?>\n" + pretty_xml_as_string + "\n"

		#print(pretty_xml_as_string)
		output.write(pretty_xml_as_string)
		output.close()



def main(argv):
	inputfile = ''
	output_dir = './'
	try:
		opts, args = getopt.getopt(argv,"ho:i:",["input=", "output-dir="])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			usage()
			sys.exit()
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
	#parseAnnotation(inputfile, output_dir, output_filepath, use_log)
	extractAnnotationFiles(inputfile, output_dir)

def usage():
	print (__doc__)
if __name__ == "__main__":
	main(sys.argv[1:])
