#CC = icc
#CC = gcc-4.0
CC = gcc
#CFLAGS = -O3 -mp -axW -tpp7 -ipo
CFLAGS = -O3 -Wall -Wextra
TARGET = stic.exe
OBJS = main.o m_alloc.o frameworks.o toolbox.o nrutil.o init.o plot.o numerical.o track.o get_field.o lin_para.o bricabrac.o othercodes_tools.o
LIBS = -lm -lClapack
#LIBS = -lm
#LIBS = -lm -lClapack -lCblas -lF77 -lg2c
#LIBS = -lm -lClapack -lCblas -lf77lapack 
#LIBS = -lm -lClapack -lCblas -lf77lapack -lg2c 
##

all :	$(TARGET)

clean :
	-rm -f $(TARGET) $(OBJS) *~ B%#* *.dat *.out *.il core

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $<

##
main.o: main.h common.h constants.h
m_alloc.o: m_alloc.h 
frameworks.o: frameworks.h 
toolbox.o: toolbox.h 
nrutil.o: nrutil.h 
init.o: init.h 
plot.o: plot.h 
numerical.o: numerical.h 
track.c: track.h 
get_field.o: get_field.h 
lin_para.o: lin_para.h 
bricabrac.o: bricabrac.h 
othercodes_tools.o: othercodes_tools.h
