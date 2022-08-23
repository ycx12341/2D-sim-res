CC            := cc
RM_UNIX       := rm -rf
RM_WIN        := del

DFLAG         := -g
OFLAG         := -o
MFLAG         := -lm
GFLAGS        := -Wall -Wextra
CFLAGS        := $(DFLAG) $(GFLAGS) -c
LFLAGS        := $(DFLAG) $(GFLAGS)

COLL_DIR      := collection
MATH_DIR      := math
MTRX_DIR      := matrix

COLL_TARGET    = arrayList \
                 collection
MATH_TARGET    = seed \
				 sample
MTRX_TARGET    = matrix

.PHONY: all coll math mtrx clean_unix

all: coll math mtrx

coll:
	for tgt in $(COLL_TARGET) ; do \
		$(CC) $(CFLAGS) $(COLL_DIR)/$$tgt.c $(OFLAG) $(COLL_DIR)/$$tgt.o $(MFLAG) ; \
	done

math:
	for tgt in $(MATH_TARGET) ; do \
		$(CC) $(CFLAGS) $(MATH_DIR)/$$tgt.c $(OFLAG) $(MATH_DIR)/$$tgt.o $(MFLAG) ; \
	done

mtrx:
	for tgt in $(MTRX_TARGET) ; do \
		$(CC) $(CFLAGS) $(MTRX_DIR)/$$tgt.c $(OFLAG) $(MTRX_DIR)/$$tgt.o $(MFLAG) ; \
	done

clean_unix:
	$(RM_UNIX) $(COLL_DIR)/*.o
	$(RM_UNIX) $(MATH_DIR)/*.o
	$(RM_UNIX) $(MTRX_DIR)/*.o

clean_win:
	$(RM_WIN) $(COLL_DIR)\*.o
	$(RM_WIN) $(MATH_DIR)\*.o
	$(RM_WIN) $(MTRX_DIR)\*.o