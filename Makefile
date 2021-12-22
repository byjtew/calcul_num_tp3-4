##########################################
# Makefile                               #
# Makefile for the code developed in TP1 #
#                                        #
# T. Dufaud                              #
##########################################
################################
# Variables for this makefile
################################
# 
CC=gcc

# 
# -- Compiler Option
#
OPTC=-O3 -Wall -fomit-frame-pointer -fPIC -mavx -DAdd_ -DF77_INTEGER=int -DStringSunStyle -fopenmp

#
# -- Directories
TPDIR=.
TPDIRSRC=$(TPDIR)/src
TPDIRTESTS=$(TPDIR)/tests

#
# -- librairies
LIBS=-llapacke -lblas -lm

# -- Include directories
INCLBLASLAPACK= -I /usr/include

INCL= -I $(TPDIR)/include $(INCLBLASLAPACK) 
#
#################################################################
# makefile
############
#
OBJENV= tp_env.o
OBJTP2ITER= lib_poisson1D.o tp2_poisson1D_iter.o
OBJTP2DIRECT= lib_poisson1D.o tp2_poisson1D_direct.o
#

all: bin/tp_testenv bin/tp2poisson1D_iter bin/tp2poisson1D_direct bin/tp4_ex5 bin/dgbmv check

testenv: bin/tp_testenv

tp2poisson1D_iter: bin/tp2poisson1D_iter

tp2poisson1D_direct: bin/tp2poisson1D_direct

tp_env.o: $(TPDIRSRC)/tp_env.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp_env.c 

dgbmv.o: $(TPDIRSRC)/dgbmv.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/dgbmv.c

lib_poisson1D.o: $(TPDIRSRC)/lib_poisson1D.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/lib_poisson1D.c 

tp2_poisson1D_iter.o: $(TPDIRSRC)/tp2_poisson1D_iter.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_iter.c  

tp2_poisson1D_direct.o: $(TPDIRSRC)/tp2_poisson1D_direct.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp2_poisson1D_direct.c  

bin/tp_testenv: $(OBJENV) 
	$(CC) -o bin/tp_testenv $(OPTC) $(OBJENV) $(LIBS)

bin/tp2poisson1D_iter: $(OBJTP2ITER)
	$(CC) -o bin/tp2poisson1D_iter $(OPTC) $(OBJTP2ITER) $(LIBS)

bin/tp2poisson1D_direct: $(OBJTP2DIRECT)
	$(CC) -o bin/tp2poisson1D_direct $(OPTC) $(OBJTP2DIRECT) $(LIBS)


bin/tp4_ex5: $(TPDIRSRC)/TP4_ex5.c
	$(CC) -o bin/tp4_ex5 $(TPDIRSRC)/TP4_ex5.c $(OPTC)

bin/tp4_ex1: $(TPDIRSRC)/TP4_ex1.sci
	scilab-cli -f $(TPDIRSRC)/TP4_ex1.sci -quit

bin/tp5_ex5: $(TPDIRSRC)/TP4_ex1.sci
	scilab-cli -f $(TPDIRSRC)/TP5_ex5.sci -quit

bin/dgbmv: dgbmv.o
	$(CC) -o bin/dgbmv $(OPTC) dgbmv.o $(LIBS)

# TESTS

check: 
	make poisson_row_col_major

poisson_row_col_major.o: $(TPDIRTESTS)/poisson_row_col_major.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRTESTS)/poisson_row_col_major.c  

poisson_row_col_major: lib_poisson1D.o poisson_row_col_major.o
	$(CC) -o bin/poisson_row_col_major $(OPTC) lib_poisson1D.o poisson_row_col_major.o $(LIBS)


# == #

run_testenv:
	bin/tp_testenv

run_tp2poisson1D_iter:
	bin/tp2poisson1D_iter

run_tp2poisson1D_direct:
	bin/tp2poisson1D_direct

clean:
	rm *.o bin/*  *.dat
