#/*************************************************************************
#*
#*                     MAKEFILE FOR SPIN.C                                
#*                                                                         
#*************************************************************************/



#######################################################################
#                                                                     #
# 1. Specify C compiler and ANSI option:                              #
#                                                                     #      
####################################################################### 

# Replace with the name of your c compiler if you don't use gcc
CC=gcc

#/*************************************************************************
#*                     SPECIFY GRID SIZE
#*************************************************************************/

#STANDARD
#HIGHSIZE=-DMDIV=65 -DSDIV=129

#HIGH
#HIGHSIZE=-DMDIV=101 -DSDIV=201

#VERY HIGH
HIGHSIZE=-DMDIV=201 -DSDIV=401


#/*************************************************************************
#*                     COMPILING FLAGS
#*************************************************************************/


# DEBUGGING OPTION
MY_OWN =-g3 -Wall 
#MY_OWN = -g3 -Wall -DDEBUG

#MY_OWN =-O3

#/*************************************************************************
#*                    SOURCE AND OBJECT MACROS
#*************************************************************************/


SOBJ=quasikerr.o findmodel.o equil.o equil_util.o nrutil.o quadrupole.o surface.o interp.o metric.o

#/*************************************************************************
#*                    MAIN COMPILING INSTRUCTIONS
#*************************************************************************/

quasik: $(SOBJ)
	$(CC) $(MY_OWN) -lm  $(HIGHSIZE) -o quasik $(SOBJ)

quasikerr.o: consts.h struct.h nrutil.h equil.h equil_util.h findmodel.h quadrupole.h surface.h interp.h metric.h quasikerr.c makefile
	$(CC) -c $(MY_OWN) $(CFLAGS) $(COPTFLAGS) $(HIGHSIZE)  quasikerr.c 

findmodel.o: consts.h struct.h nrutil.h equil.h equil_util.h surface.h findmodel.h findmodel.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(HIGHSIZE)   findmodel.c

equil.o:equil.h equil_util.h nrutil.h consts.h equil.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(HIGHSIZE)   equil.c

equil_util.o:equil_util.h nrutil.h consts.h equil_util.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(HIGHSIZE)   equil_util.c

nrutil.o:nrutil.h nrutil.c makefile
	$(CC) -c $(COPTFLAGS)   nrutil.c

surface.o: consts.h struct.h nrutil.h equil.h equil_util.h findmodel.h surface.h quadrupole.h surface.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(HIGHSIZE)  surface.c 

quadrupole.o: consts.h struct.h nrutil.h equil.h equil_util.h findmodel.h surface.h interp.h quadrupole.h quadrupole.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(HIGHSIZE)  quadrupole.c

interp.o: consts.h nrutil.h equil_util.h interp.h interp.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS) $(HIGHSIZE)  interp.c

metric.o: consts.h metric.c makefile
	$(CC) -c $(MY_OWN) $(COPTFLAGS)  metric.c
