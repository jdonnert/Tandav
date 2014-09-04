pro softening

	common globals, tandav, cosmo

	N = 1024
	rmin = 1e-4
	rmax = 1d4

	i = lindgen(N)
	
	dr = alog10(rmax/rmin) / (N-1)
	r = rmin * 10D^(i * dr)

	eps = 1D

	h = 3/eps; 

	kernel = make_array(N, /double, val=0)
	force = make_array(N, /double, val=0)
	poten = make_array(N, /double, val=0)
	
	for i = 0,N-1 do begin

		kernel[i] = WC2(r[i], h)

		force[i] = force_WC2(1,1,r[i], eps)
		
		poten[i] = potential_WC2(1, r[i], eps)
	end

	!p.multi[1] = 3

	plot, r, kernel, xrange=[rmin, 1.1*h], xstyle=1, xtitle='r [cm]', $
		ytitle=textoidl('\rho [g/cm^3]')
	plot, r, force, /ylog, /xlog, xtitle='r [cm]', ytitle='F [g cm/s^2]'
	plot, r, -poten, /ylog, /xlog, xtitle='r [cm]', $
		ytitle=textoidl('\phi [g cm^2/s^2]')

	stop

	return 
end

function WC2, r, h

	u = r/h

	if u gt 1 then $
		return, 0

	return, 21D/(2D*!pi) /h^3 * (1-u)^4 * (1+4*u)
end

function force_WC2, M0, M1, r, eps

	grav_const = 6.6720000e-08

	h = 3/eps; 

	u = r/h

	if r gt h then $
		u = 1 

	force_law = (14*u^3 - 84*u^5 + 140*u^6 - 90*u^7 + 21*u^8)/ r^2

	return, grav_const * M0 * M1  * force_law
end

function potential_WC2, M, r ,eps

	grav_const = 6.6720000e-08

	h = 3/eps 

	if r gt h then $
		return, -grav_const * M  / r

	u = double(r/h)

	wc_poly = (7*u^3 - 21*u^5 + 28*u^6 - 15*u^7 + 3*u^8 - 3D*u)

	if u le 1 then print, u, r, h,wc_poly

	return,  grav_const * M / r * wc_poly
end
