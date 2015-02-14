pro make_ICs

	common  globals, tandav, cosmo

	msol = 1.989d33

	seed = 14041981L

	npart = 100000L

	mass = 1d15 * Msol / tandav.mass ; code units 
	a_hernq = 924D 		

	mpart = mass / npart

	; positions

	pos = make_array(3, npart, /double, val=0)
	
	sqrt_q = sqrt(randomu(seed, npart)) * 0.97

	r = a_hernq  * sqrt_q / (1-sqrt_q) 

	theta = acos(2 * randomu(seed, npart) - 1)
	phi = 2*!pi * randomu(seed, npart) 

	pos[0,*] = r * sin(theta) * cos(phi)
	pos[1,*] = r * sin(theta) * sin(phi)
	pos[2,*] = r * cos(theta) 

	print, 'npart = '+strn(npart), minmax(r) 
	print, 'bounds', minmax(pos[0,*]), minmax(pos[1,*]), minmax(pos[2,*])

	plot, pos[0,*], pos[1,*], /iso, psym=3

	; velocities

	vel = make_array(3, npart, /double, val=0)

	for i = 0, npart-1 do begin

		r = sqrt(pos[0,i]^2 + pos[1,i]^2 + pos[2,i]^2)

		pot = hernquist_potential(r, a_hernq, mass)

 		vmax = sqrt(-2*pot)
 		Emax = -pot
 		qmax = 4*!pi*vmax^2/mass $
 			* hernquist_distribution_function(Emax, a_hernq, mass)

 		while 1 do begin ; rejection sampling
	
 			lower_bound = qmax * randomu(seed)

 			v = vmax * randomu(seed)

 			Etot = 0.5*v^2 + pot

 			q = 4*!pi * v^2/mass $
 				* hernquist_distribution_function(-Etot, a_hernq, mass)

 			if q ge lower_bound then $
 				break
 		end
 
 		theta = acos(2 * randomu(seed) - 1)
 		phi =  2 * !pi * randomu(seed)

 		vel[0,i] = v * sin(theta) * cos(phi)
 		vel[1,i] = v * sin(theta) * sin(phi)
 		vel[2,i] = v * cos(theta)
 	end
	
	; output

	head = tandav.MakeHead()

	head.npart[1] = LONG(npart)
	head.parttotal = head.npart
	head.massarr[1] = mpart
	head.time = 0
	head.redshift = 0
	head.num_files = 1
	head.boxsize = 1d15
	head.Omega0 = 1D
	head.OmegaLambda = 0.7D
	head.HubbleParam = 0.7200D

	fname = './IC_Hernquist_Halo'

	print, 'IC File', fname

	tandav.WriteHead, fname, head
	
	tandav.AddBlock, fname, float(pos), 'POS'
	tandav.AddBlock, fname, float(vel), 'VEL'
	tandav.AddBlock, fname, ulong(indgen(npart)+1), 'ID'

	return
end

function hernquist_potential, r, a, mass

	common  globals, tandav, cosmo

	return, -tandav.Grav * mass / (r+a)
end

function hernquist_distribution_function, E, a, mass

	common  globals, tandav, cosmo
	
	prefac = 1D / (sqrt(2) * (2*!pi)^3 * (tandav.Grav*mass*a)^(1.5) )

	q2 = a * E / (tandav.Grav * mass)

	fE = prefac * mass * sqrt(q2) / (1-q2)^2 $
		*( (1 - 2*q2) * (8*q2*q2 - 8*q2 - 3) $
		+ (3*asin(sqrt(q2)))/sqrt(q2*(1-q2)) )

	return, fE
end
