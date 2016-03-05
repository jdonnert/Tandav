pro show_box, snap

	common globals, gadget, tandav, cosmo

	if not keyword_set(snap) then $
		snap = 000

	fname = 'snap_'+strn(snap, len=3, padc='0')

	pos = tandav.readsnap(fname, 'POS', head=head)

	pos[0,*] += head.boxsize/2
	pos[1,*] += head.boxsize/2
	pos[2,*] += head.boxsize/2

	bad = where(pos[0,*] gt head.boxsize)
	pos[0,bad] -= head.boxsize
	bad = where(pos[1,*] gt head.boxsize)
	pos[1,bad] -= head.boxsize
	bad = where(pos[2,*] gt head.boxsize)
	pos[2,bad] -= head.boxsize
	
	good = where(pos[2,*] lt 1d40 )

	plot, pos[0,good], pos[1,good], psym=3, /iso, xtitle="x [kpc]", $
		ytitle="y [kpc]"

	return
end


pro make_frames

	common globals, tandav, cosmo

	cosmo.set, 1

	outpath = "./frames/"

	nsnap = 100L

	win, x=1024, y=1024

	for snap = 0, nsnap-1 do begin

		show_box, snap

		outname = outpath+'frame_'+strn(snap, len=3, padc='0')

		save_screen, outname

	end
	
	return
end
