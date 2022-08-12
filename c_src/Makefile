CC            := cc
RM_UNIX       := rm -rf

DFLAG         := -g
OFLAG         := -o
MFLAG         := -lm
GFLAGS        := -Wall -Wextra
CFLAGS        := $(DFLAG) $(GFLAGS) -c
LFLAGS        := $(DFLAG) $(GFLAGS)

ALGO_DIR      := algorithm
UTIL_DIR      := util
UTIL_COLL_DIR := $(UTIL_DIR)/collection
UTIL_MATH_DIR := $(UTIL_DIR)/math
UTIL_MTRX_DIR := $(UTIL_DIR)/matrix

EXE           := main

.PHONY: all algorithm util clean_unix main

all: clean_unix algorithm util main

algorithm:
	cd $(ALGO_DIR) && make

util:
	cd $(UTIL_DIR) && make

main:
	$(CC) $(LFLAGS)          \
	    $(ALGO_DIR)/*.o      \
		$(UTIL_COLL_DIR)/*.o \
		$(UTIL_MATH_DIR)/*.o \
		$(UTIL_MTRX_DIR)/*.o \
		$(EXE).c             \
		$(OFLAG) $(EXE) $(MFLAG)

clean_unix:
	$(RM_UNIX) $(EXE)
	cd $(ALGO_DIR) && make clean_unix
	cd $(UTIL_DIR) && make clean_unix