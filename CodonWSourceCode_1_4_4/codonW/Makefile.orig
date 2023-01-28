override cflags  = $(CFLAGS) -g

objects  = codon_us.o codons.o open_fil.o commline.o menu.o tester.o coresp.o
linked   = rscu cu aau raau tidy reader cutab cutot transl bases base3s dinuc cai fop gc3s gc cbi enc

CC=cc
CFLAGS= -O -DBSD
LN=ln -f


all: codonw links   

codonw: $(objects)
	$(CC) $(CFLAGS)  $(objects) -o codonw -lm

clean:
	\rm -f $(objects)

cleanall:
	\rm -f $(objects) codonw Makefile $(linked)

realclean:
	\rm -f $(objects) codonw Makefile $(linked)

codon_us.o: codon_us.c codonW.h 
	$(CC) -c $(CFLAGS) codon_us.c  

menu.o: menu.c codonW.h 
	$(CC) -c $(CFLAGS) menu.c

codons.o: codons.c codonW.h 
	$(CC) -c $(CFLAGS) codons.c

coresp.o: coresp.c codonW.h 
	$(CC) -c $(CFLAGS) coresp.c

open_fil.o:    open_fil.c codonW.h
	$(CC) -c $(CFLAGS) open_fil.c

commline.o:    commline.c codonW.h 
	$(CC) -c $(CFLAGS) commline.c

tester.o:      tester.c codonW.h
	$(CC) -c $(CFLAGS) tester.c

links: codonw
		$(LN) codonw rscu
		$(LN) codonw cu
		$(LN) codonw aau
		$(LN) codonw raau
		$(LN) codonw tidy
		$(LN) codonw reader
		$(LN) codonw cutab
		$(LN) codonw cutot
		$(LN) codonw transl
		$(LN) codonw bases
		$(LN) codonw base3s
		$(LN) codonw dinuc
		$(LN) codonw cai
		$(LN) codonw fop
		$(LN) codonw gc3s
		$(LN) codonw gc
		$(LN) codonw cbi
		$(LN) codonw enc

