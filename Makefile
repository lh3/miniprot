CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O3
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		kalloc.o kthread.o sys.o misc.o table.o options.o ntseq.o index.o
PROG=		miniprot
LIBS=		-lz

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

miniprot:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) fmd-occ *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

index.o: mppriv.h miniprot.h kalloc.h kvec-km.h kthread.h
kalloc.o: kalloc.h
kthread.o: kthread.h
main.o: miniprot.h ketopt.h
misc.o: mppriv.h miniprot.h ksort.h
ntseq.o: mppriv.h miniprot.h kalloc.h kseq.h
options.o: miniprot.h
sys.o: mppriv.h miniprot.h
table.o: miniprot.h
