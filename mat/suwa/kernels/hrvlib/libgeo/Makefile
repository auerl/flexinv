# moved this library to hrvsrc and added some new flags for the new
# compiler (-O3 -cg89 -libmil).  g.e. 2/91
#
LIB = libgeo
FFLAGS = $(FFLAG)

.f.a:
	f77 -c $(FFLAGS) $<
	ar rv $@ $*.o
	rm -f $*.o

.c.a:
	cc -c -O $<
	ar rv $@ $*.o
	rm -f $*.o

$(HRVLIB)/$(LIB).a: \
	$(HRVLIB)/$(LIB).a(ahouse.o) \
	$(HRVLIB)/$(LIB).a(asrins.o) \
	$(HRVLIB)/$(LIB).a(sroins.o) \
	$(HRVLIB)/$(LIB).a(bhouse.o) \
	$(HRVLIB)/$(LIB).a(bitrv.o) \
	$(HRVLIB)/$(LIB).a(convsp.o) \
	$(HRVLIB)/$(LIB).a(cool2.o) \
	$(HRVLIB)/$(LIB).a(ddabs.o) \
	$(HRVLIB)/$(LIB).a(delaz.o) \
	$(HRVLIB)/$(LIB).a(dfloat.o) \
	$(HRVLIB)/$(LIB).a(dot.o) \
	$(HRVLIB)/$(LIB).a(drsple.o) \
	$(HRVLIB)/$(LIB).a(drspln.o) \
	$(HRVLIB)/$(LIB).a(dspint.o) \
	$(HRVLIB)/$(LIB).a(fixrl.o) \
	$(HRVLIB)/$(LIB).a(float2.o) \
	$(HRVLIB)/$(LIB).a(four2.o) \
	$(HRVLIB)/$(LIB).a(instr.o) \
	$(HRVLIB)/$(LIB).a(legndr.o) \
	$(HRVLIB)/$(LIB).a(matinv.o) \
	$(HRVLIB)/$(LIB).a(rfour.o) \
	$(HRVLIB)/$(LIB).a(rsple.o) \
	$(HRVLIB)/$(LIB).a(rspln.o) \
	$(HRVLIB)/$(LIB).a(spint.o) \
	$(HRVLIB)/$(LIB).a(tfr.o) \
	$(HRVLIB)/$(LIB).a(tfg.o) \
	$(HRVLIB)/$(LIB).a(wier.o) \
	$(HRVLIB)/$(LIB).a(trnslt.o) \
	$(HRVLIB)/$(LIB).a(remavl.o) \
	$(HRVLIB)/$(LIB).a(tapclr.o)

$(HRVLIB)/$(LIB).a: ; ranlib $(HRVLIB)/$(LIB).a
