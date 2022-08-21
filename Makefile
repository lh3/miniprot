CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O3
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		kalloc.o kthread.o sys.o misc.o table.o options.o ntseq.o sketch.o \
			index.o bseq.o chain.o nasw-s.o hit.o format.o map.o
PROG=		miniprot
LIBS=		-lpthread -lz

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl -lm
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

bseq.o: kvec-km.h kalloc.h mppriv.h miniprot.h kseq.h
chain.o: mppriv.h miniprot.h kseq.h kalloc.h
format.o: kseq.h mppriv.h miniprot.h
hit.o: mppriv.h miniprot.h kseq.h kalloc.h
index.o: mppriv.h miniprot.h kseq.h kalloc.h kvec-km.h kthread.h
kalloc.o: kalloc.h
kthread.o: kthread.h
main.o: mppriv.h miniprot.h kseq.h ketopt.h
map.o: mppriv.h miniprot.h kseq.h kthread.h kvec-km.h kalloc.h
misc.o: mppriv.h miniprot.h kseq.h ksort.h
nasw-s.o: nasw.h kalloc.h mppriv.h miniprot.h kseq.h
ntseq.o: mppriv.h miniprot.h kseq.h kalloc.h
options.o: miniprot.h
sketch.o: mppriv.h miniprot.h kseq.h kalloc.h kvec-km.h
sys.o: mppriv.h miniprot.h kseq.h
table.o: miniprot.h
