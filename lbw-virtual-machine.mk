#######################################
# lbw-virtual-machine.mk
# Default options for lbw-virtual-machine computer
#######################################
CC=gcc
LIBSLOCAL=-L/usr/lib -llapack -lblas -lm
INCLUDEBLASLOCAL=-I/usr/include
OPTCLOCAL=-fPIC -march=native
