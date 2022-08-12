CC            := cc
RM_UNIX       := rm -rf

DFLAG         := -g
OFLAG         := -o
MFLAG         := -lm
GFLAGS        := -Wall -Wextra
CFLAGS        := $(DFLAG) $(GFLAGS) -c
LFLAGS        := $(DFLAG) $(GFLAGS)

TARGET 		  := calculate_scc \
                 pars

.PHONY: all clean_unix

all:
	for tgt in $(TARGET) ; do \
		$(CC) $(CFLAGS) $$tgt.c $(OFLAG) $$tgt.o $(MFLAG) $(MFLAG) ; \
	done

clean_unix:
	$(RM_UNIX) *.o