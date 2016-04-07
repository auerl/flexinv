	subroutine print_bldcommand()
c  print options for the build program
	print*, 'Usage:'
	print*, '--------------------------------------------------'
	print*, 'build_ata_zetavbar [opt]'
	print*, 'opt:	(options in constructing ATA)'
	print*, '1  ----  normal build for SV/PV & SH/PH ATA'
	print*, '2  ----  build Vdiff and V_voigt from V_xh & V_xv'
	print*, '	  xh and xv could be P or S'
	print*, '	  Vdiff=Vxv-Vxh and V_avg=(Vxh+2Vxv)/3'
	print*, '3  ----  build Vdiff and V_avg from V_xh & V_xv'
	print*, '	  xh and xv could be P or S'
	print*, '	  Vdiff=Vxv-Vxh and V_avg=Vxh+Vxv'
	print*, '--------------------------------------------------'
	stop
	return
	end
