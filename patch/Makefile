CC=			g++
CFLAGS=		-g -Wall -O2 -Wno-all #-Wc++-compat 
CPPFLAGS=	-Wno-write-strings -Wno-deprecated-declarations
INCLUDES=	-I.
OBJS=		sys.o sdict.o paf.o asg.o common.o hit.o asm.o rsig.o rutils.o kalloc.o rawasm.o main.o 
PROG=		miniasm minidot
LIBS=		-lm -lz -lpthread#-luuid

# Rawasm-related directories
RA_LIB_DIR = ./rawasm/lib
RA_INCLUDE_DIR = ./rawasm/include
RA_SRC_DIR = ./rawasm/src
INCLUDES+=-I${RA_INCLUDE_DIR}

# HDF5 Libraries & Headers
HDF5_INCLUDE_DIR ?= ${RA_INCLUDE_DIR}/hd5
INCLUDES+=-I${HDF5_INCLUDE_DIR}
LIBS+=${RA_LIB_DIR}/libhdf5.a 

# POD5 Libraries & Headers
POD5_INCLUDE_DIR = ${RA_INCLUDE_DIR}/pod5
INCLUDES+=-I${POD5_INCLUDE_DIR}
LIBS+=${RA_LIB_DIR}/libpod5_format.a ${RA_LIB_DIR}/libarrow.a ${RA_LIB_DIR}/libjemalloc_pic.a ${RA_LIB_DIR}/libzstd.a

# SLOW5 Libraries & Headers
SLOW5_INCLUDE_DIR ?= ${RA_INCLUDE_DIR}/slow5
INCLUDES+=-I${SLOW5_INCLUDE_DIR}
LIBS+=${RA_LIB_DIR}/libslow5.a

# UUID Libraries & Headers
UUID_INCLUDE_DIR ?= ${RA_INCLUDE_DIR}/libuuid
INCLUDES+=-I${UUID_INCLUDE_DIR}
LIBS+=${RA_LIB_DIR}/libuuid.a


all:$(PROG)

miniasm_obj:
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) sys.c -o sys.o
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) asg.c -o asg.o
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) asm.c -o asm.o
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) common.c -o common.o
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) hit.c -o hit.o
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) dotter.c -o dotter.o
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) paf.c -o paf.o
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) sdict.c -o sdict.o


rawasm_obj:
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $(RA_SRC_DIR)/rsig.c -o rsig.o
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $(RA_SRC_DIR)/rutils.c -o rutils.o
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $(RA_SRC_DIR)/kalloc.c -o kalloc.o
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $(RA_SRC_DIR)/rawasm.c -o rawasm.o

main_obj:
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) main.c -o main.o


miniasm:miniasm_obj rawasm_obj main_obj
		$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS)

minidot:paf.o sdict.o dotter.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)
