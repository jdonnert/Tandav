; particle on Kepler orbit tests time integration scheme
;  Binney & Tremaine pp148
pro make_ICs
	
	common globals, gadget, cosmo

	M = [1d11, 1d0] / gadget.mass 

	e = 0.9 ; eccentricity

	x0 = 200D

	eta = 0 ; eccentric anomaly
	psi0 = -2 * atan( sqrt( (1+e) / (1-e) ) * tan(0.5*eta)  ) ; eq. 3.29

	pos = [ [0, 0, 0], [x0, 0, 0] ]

	r = sqrt( pos[0,1]^2 + pos[1,1]^2 + pos[2,1]^2 )
	
	v = sqrt( r * gadget.grav * M[0] * (1 - e) / (1-e*cos(eta)) ) ; Keplers eq.
	vel = [ [0,0,0], [0,0,-v] ] ; L = r x v in y direction

	L2 = (pos[0,1]*vel[2,1])^2 ; square of angular momentum L = r x vel

	C  = e/L2 * gadget.grav*M[0]
	a = L2 / (gadget.grav*M[0]*(1-e^2))

	nstep = 1000
	x = make_array(nstep)
	y = make_array(nstep)

	for i=0,nstep-1 do begin

		psi = i*2*!pi / nstep
	
		r = a*(1-e^2) / (1 + e*cos(psi-psi0))

		x[i] = r * cos(psi)
		y[i] = r * sin(psi)

	end

	plot, x, y, /iso
	oplot, x[0]*[1,1], y[0]*[1,1], psym=4
	oplot, [0,0], [0,0], psym=1

	fout = '../IC_Kepler'

	print, fout

	head = gadget.MakeHead()
	head.npart = [2, 0, 0, 0, 0, 0]
	head.massarr[*] = 0
	head.time = 0
	head.redshift = 0
	head.parttotal = head.npart

	gadget.WriteHead, fout, head
	gadget.AddBlock, fout, float(pos), 'POS'
	gadget.AddBlock, fout, float(vel), 'VEL'
	gadget.AddBlock, fout, ulong([1,2]), 'ID'
	gadget.AddBlock, fout, float(M), 'MASS'

	return
end
