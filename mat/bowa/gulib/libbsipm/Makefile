FFLAGS= -extend_source
#FFLAGS=  -e -O3 -cg89 -libmil 
#CFLAGS = -O -cg89 -libmil -Bstatic
CFLAGS = -O -libmil -Bstatic
#CFLAGS= -g
#FFLAGS= -g

LIBDIR =../
LIB = libbsipm

.f.a:
	ifort -c $(FFLAGS) $<
	ar rv $@ $*.o
	rm -f $*.o

.c.a:
	cc -c $(CFLAGS)  $<
	ar rv $@ $*.o
	rm -f $*.o

$(LIBDIR)/$(LIB).a: \
	$(LIBDIR)/$(LIB).a(addarow.o) \
	$(LIBDIR)/$(LIB).a(add_dcshift.o) \
	$(LIBDIR)/$(LIB).a(addmat_spline.o) \
	$(LIBDIR)/$(LIB).a(addatax_lower.o) \
	$(LIBDIR)/$(LIB).a(conv2geocen.o) \
	$(LIBDIR)/$(LIB).a(dampdiag.o) \
	$(LIBDIR)/$(LIB).a(damp_norm_grad.o) \
	$(LIBDIR)/$(LIB).a(dampdiag_topo.o) \
	$(LIBDIR)/$(LIB).a(dot.o) \
	$(LIBDIR)/$(LIB).a(dot_splines.o) \
	$(LIBDIR)/$(LIB).a(filesubs.o) \
	$(LIBDIR)/$(LIB).a(find_rsplit.o) \
	$(LIBDIR)/$(LIB).a(getataf_lower.o) \
	$(LIBDIR)/$(LIB).a(get_dampfac.o) \
	$(LIBDIR)/$(LIB).a(handle_2mod.o) \
	$(LIBDIR)/$(LIB).a(initmat.o) \
	$(LIBDIR)/$(LIB).a(modify_arow.o) \
	$(LIBDIR)/$(LIB).a(neweuler.o) \
	$(LIBDIR)/$(LIB).a(print_bldcommand.o) \
	$(LIBDIR)/$(LIB).a(print_ipmavg.o) \
	$(LIBDIR)/$(LIB).a(print_ipmstat.o) \
	$(LIBDIR)/$(LIB).a(readdmat2d.o) \
	$(LIBDIR)/$(LIB).a(readheader_arow.o) \
	$(LIBDIR)/$(LIB).a(readheader.o) \
	$(LIBDIR)/$(LIB).a(read_spline.o) \
	$(LIBDIR)/$(LIB).a(read_spline1.o) \
	$(LIBDIR)/$(LIB).a(read_topo.o) \
	$(LIBDIR)/$(LIB).a(readsumres.o) \
	$(LIBDIR)/$(LIB).a(readata.o) \
	$(LIBDIR)/$(LIB).a(stat_lower.o) \
	$(LIBDIR)/$(LIB).a(writeata.o) \
	$(LIBDIR)/$(LIB).a(write_bsmodel.o) \
	$(LIBDIR)/$(LIB).a(writeheader.o) \
	$(LIBDIR)/$(LIB).a(writeipm_c.o)

$(LIBDIR)/$(LIB).a: ; ranlib $(LIBDIR)/$(LIB).a
