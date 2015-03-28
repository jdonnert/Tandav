pro show_box, snap

	common globals, tandav, cosmo

	cosmo.set, 1

	if not keyword_set(snap) then $
		snap = 000

	fname = 'snap_'+strn(snap, len=3, padc='0')

	pos = tandav.readsnap(fname, 'POS', head=head)
	
	plot, pos[0,*], pos[1,*], psym=3, /iso

	good = where(pos[2,*] lt 1e4 )

	plot, pos[0,good], pos[1,good], psym=3, /iso

	return
end
