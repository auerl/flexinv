FFLAGS= -extend_source  
CFLAGS = -O -cg89 -libmil -Bstatic

LIBDIR = ../
LIB = libkernel

.f.a:
	ifort -c $(FFLAGS) $<
	ar rv $@ $*.o
	rm -f $*.o

.c.a:
	cc -c $(CFLAGS)  $<
	ar rv $@ $*.o
	rm -f $*.o

$(LIBDIR)/$(LIB).a: \
	$(LIBDIR)/$(LIB).a(addpath_kernel.o) \
	$(LIBDIR)/$(LIB).a(addray.o) \
	$(LIBDIR)/$(LIB).a(addsurfspline.o) \
	$(LIBDIR)/$(LIB).a(anisoker_ttpath.o) \
	$(LIBDIR)/$(LIB).a(anisoker_ttpath_topo.o) \
	$(LIBDIR)/$(LIB).a(bsplinefun.o) \
	$(LIBDIR)/$(LIB).a(calcdrdl.o) \
	$(LIBDIR)/$(LIB).a(calc_fwdtransmat.o) \
	$(LIBDIR)/$(LIB).a(calcqvec.o) \
	$(LIBDIR)/$(LIB).a(calc_tessa.o) \
	$(LIBDIR)/$(LIB).a(cder.o) \
	$(LIBDIR)/$(LIB).a(check_bot.o) \
	$(LIBDIR)/$(LIB).a(circle.o) \
	$(LIBDIR)/$(LIB).a(conv2geocen.o) \
	$(LIBDIR)/$(LIB).a(deriv.o) \
	$(LIBDIR)/$(LIB).a(diffract_path.o) \
	$(LIBDIR)/$(LIB).a(dist.o) \
	$(LIBDIR)/$(LIB).a(dist_phi.o) \
	$(LIBDIR)/$(LIB).a(drdzeta.o) \
	$(LIBDIR)/$(LIB).a(drspln.o) \
	$(LIBDIR)/$(LIB).a(drsple.o) \
	$(LIBDIR)/$(LIB).a(dxt.o) \
	$(LIBDIR)/$(LIB).a(f.o) \
	$(LIBDIR)/$(LIB).a(fdelta.o) \
	$(LIBDIR)/$(LIB).a(find_impactpt.o) \
	$(LIBDIR)/$(LIB).a(fqanis.o) \
	$(LIBDIR)/$(LIB).a(fqanis_kernel.o) \
	$(LIBDIR)/$(LIB).a(gauslv.o) \
	$(LIBDIR)/$(LIB).a(getbsreg.o) \
	$(LIBDIR)/$(LIB).a(getseg_path.o) \
	$(LIBDIR)/$(LIB).a(getkIn_int.o) \
	$(LIBDIR)/$(LIB).a(getmod.o) \
	$(LIBDIR)/$(LIB).a(getmod_iso.o) \
	$(LIBDIR)/$(LIB).a(getq.o) \
	$(LIBDIR)/$(LIB).a(getpart.o) \
	$(LIBDIR)/$(LIB).a(ichdec.o) \
	$(LIBDIR)/$(LIB).a(ifbetween.o) \
	$(LIBDIR)/$(LIB).a(ifsecondary.o) \
	$(LIBDIR)/$(LIB).a(integr.o) \
	$(LIBDIR)/$(LIB).a(openfil_easy.o) \
	$(LIBDIR)/$(LIB).a(prray_kernel.o) \
	$(LIBDIR)/$(LIB).a(qtau.o) \
	$(LIBDIR)/$(LIB).a(radcontribute.o) \
	$(LIBDIR)/$(LIB).a(radcontribute_drdl.o) \
	$(LIBDIR)/$(LIB).a(radfun.o) \
	$(LIBDIR)/$(LIB).a(radfun_u7l5.o) \
	$(LIBDIR)/$(LIB).a(raynam.o) \
	$(LIBDIR)/$(LIB).a(raynam_topo.o) \
	$(LIBDIR)/$(LIB).a(raypath_kernel.o) \
	$(LIBDIR)/$(LIB).a(raypath_kernel_topo.o) \
	$(LIBDIR)/$(LIB).a(readmod.o) \
	$(LIBDIR)/$(LIB).a(readrad.o) \
	$(LIBDIR)/$(LIB).a(rot_euler.o) \
	$(LIBDIR)/$(LIB).a(savearow.o) \
	$(LIBDIR)/$(LIB).a(savearow_topo.o) \
	$(LIBDIR)/$(LIB).a(saveheader.o) \
	$(LIBDIR)/$(LIB).a(savepartvel.o) \
	$(LIBDIR)/$(LIB).a(spline_int.o) \
	$(LIBDIR)/$(LIB).a(surfreflker.o) \
	$(LIBDIR)/$(LIB).a(surftransker.o) \
	$(LIBDIR)/$(LIB).a(tder.o) \
	$(LIBDIR)/$(LIB).a(transform.o) \
	$(LIBDIR)/$(LIB).a(transinv.o) \
	$(LIBDIR)/$(LIB).a(writekin.o) \
	$(LIBDIR)/$(LIB).a(xlinearint.o) \
	$(LIBDIR)/$(LIB).a(xmodelvel.o) \
	$(LIBDIR)/$(LIB).a(zero.o) \
	$(LIBDIR)/$(LIB).a(zero1.o) \
	$(LIBDIR)/$(LIB).a(radcontribute_drdl_vox.o) \
	$(LIBDIR)/$(LIB).a(voxel_int.o) \
	$(LIBDIR)/$(LIB).a(gceq.o) \
	$(LIBDIR)/$(LIB).a(span.o) \
	$(LIBDIR)/$(LIB).a(isqre.o) \
	$(LIBDIR)/$(LIB).a(linint.o) \
	$(LIBDIR)/$(LIB).a(paramjsbw.o) 

$(LIBDIR)/$(LIB).a: ; ranlib $(LIBDIR)/$(LIB).a
