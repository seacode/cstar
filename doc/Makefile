# $Id: Makefile 599 2012-09-03 23:33:29Z jsibert $
#
# Author: James Ianelli
#
.PHONY: default clean

#REV=`svnversion -n ../..`
#VER=`cat ../../VERSION`


default: Doxyfile ./dox-share/mainpage.txt 
	mkdir -vp dox/html
#	Insert some linux-specific deltails into configuration file
#	See Doxygen FAQ "Can I configure doxygen from the command line?"
	(cat Doxyfile ; echo "PROJECT_NUMBER = "$(VER)"."$(REV); echo "FILE_VERSION_FILTER  = ./getID.sh"; ) | doxygen -

clean:
	rm -rf dox
