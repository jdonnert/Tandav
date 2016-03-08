; make resolution independent DM only cosmological initial conditions
; Power spectrum, Eisenstein & Hu 1998/99 
; Algorithm from N-GenIC, Springel 2006
; Zeldovich approximation, Bertschinger 1998 Ann. Rev., Efstathiou 1985
; [gadget]: use gadget units in velocity and hbpar=1
; [showPk]: plot P(k) from sampled data and analytic formula

pro make_ICs, N, boxsize=boxsize, gadget=gadget, showPk=showPk

	common parameters, tandav,  hbpar,  Omega_M, Omega_L, Tcmb

	tandav = obj_new('TANDAVCODEOBJECT')
	
	Mpc2cm = 3.0856802d24 

	if not keyword_set(N) then $
		N = 64UL $
	else $
		N = ulong64(N)

	npart = N^3

	if not keyword_set(boxsize) then $
		boxsize = 150000D $
	else $
		boxsize = double(boxsize)

	prim_idx = 1D  ; index of primordial power spectrum
	sigma8 = 0.8D  ; normalisation of P(k) at 8 Mpc
	Tcmb = 2.728D  ; temperature of the CMB [K]
	hbpar = 0.7D   ; Hubble parameter
	Omega_L = 0.7D ; Dark Energy at z=0
	Omega_M = 0.3D ; Total Matter, EH98 call this Omega_0
	Omega_B = 0.04 ; Baryons 
	
	H0 = 100D * 1d5 / Mpc2cm  * tandav.time ; tandav Units

	if not keyword_set(gadget) then $ ; hbpar = 1
		H0 *= hbpar

	set_Pk_vars, Omega_M, Omega_B, hbpar, sigma8, prim_idx, tandav.length

	pos = make_mesh(N, Boxsize)

	a_init = 1/(1 + 63D); find_initial_time(N, boxsize)
	z_init = 1D + 1D/a_init

	Da = Growth_Function(1D) / Growth_Function(a_init) ; growth factor

	displ_grid = displacement_fields(N, Boxsize, Da, showPk=showPk)

	sdev = [stddev(displ_grid[0,*,*,*]), stddev(displ_grid[1,*,*,*]), $
			stddev(displ_grid[2,*,*,*]) ]

	cellsize = boxsize/N

	displ = make_array(3, npart, val=0D)
	displ[0,*] = CIC(N, pos/cellsize, displ_grid[0,*,*,*])
	displ[1,*] = CIC(N, pos/cellsize, displ_grid[1,*,*,*])
	displ[2,*] = CIC(N, pos/cellsize, displ_grid[2,*,*,*])

	pos += displ
	
	sdev = [stddev(displ[0,*]), stddev(displ[1,*]), stddev(displ[2,*]) ]

	f0 = ( Omega_M / a_init^3D / g(a_init)^2D)^(0.6D)  ; EH98, eq. 28
	displ2vel = a_init * g(a_init) * H0 * f0; EH99, eq. 29
	
	vel2comov = 1D ; tandav comoving vel ~ 1/a^2

	if keyword_set(gadget) then $
		vel2comov = 1D/sqrt(a_init)  ; Gadget comoving vel ~ 1/sqrt(a)

	vel = displ * displ2vel * vel2comov

	constrain_particles_to_box, pos, boxsize

	id = ulindgen(npart)+1
	
	rho_crit = 3D * H0^2D/(8D*!pi*tandav.Grav) ; comoving tandav units
	total_mass = Omega_M * rho_crit * boxsize^3D

	print, "npart      = "+strn(npart)
	print, "H(0)       = "+strn(H0)+" code units"
	print, "H(a_i)     = "+strn(H0*g(a_init))
	print, "hbpar      = "+strn(hbpar)
	print, "Initial a  = "+strn(a_init)
	print, "Initial z  = "+strn(z_init)
	print, "Dplus      = "+strn(Da) ; growth factor	
	print, "F0(a_ini)  = "+strn(f0) 
	print, "displ2vel  = "+strn(displ2vel*vel2comov)
	print, "max(vel)   = "+strn(max(vel))
	print, "rho_c(z0)  = "+strn(rho_crit)
	print, "Mtot       = "+strn(total_mass)
	print, "Mpart      = "+strn(total_mass/npart)

	print, "Displacements: "
	print, "   max     = "+strn(max(displ))+" kpc/h" 
	print, "   mean    = ", mean(displ[0,*]),mean(displ[1,*]),mean(displ[2,*])
	print, "   std_dev = ", sdev 


	head = tandav.make_head()
	head.npart = [0,1,0,0,0,0] * npart
	head.massarr = [0,1,0,0,0,0] * total_mass/npart
	head.time = a_init
	head.redshift = z_init
	head.flag_sfr = 0
	head.flag_feedback = 0
	head.parttotal = head.npart
	head.flag_cooling = 0
	head.num_files = 1
	head.boxsize = boxsize
	head.omega0 = Omega_M
	head.omegalambda = Omega_L
	head.hubbleparam = hbpar

	if not keyword_set(fout) then $
		fout = "IC_Cosmo_Box_IDL_"+strn(N)

	tandav.write_head, fout, head
	tandav.add_block, fout, float(pos), "POS"
	tandav.add_block, fout, float(vel), "VEL"
	tandav.add_block, fout, ulong(id), "ID"

	epsilon = (Boxsize^3D/npart)^(1D/3D) / 7D

	min_boxsize = find_min_boxsize(N)

	print, "-------------------------------------"
	print, "fout            = ", strn(fout)
	print, "boxsize         = ", strn(boxsize)
	print, "a_init          = "+strn(a_init)
	print, "z_init          = "+strn(1/a_init-1)
	print, "Grav Softening  = "+strn(epsilon)
	print, "linear mode (sigma_box=1) is at "+strn(min_boxsize)+" unit length "

	return
end

function make_mesh, N,  boxsize
	
	pos = make_array(3, N^3, val=0D)

	fac = boxsize/N

	for i = 0UL, N-1 do begin
		for j = 0UL, N-1 do begin
			for k = 0UL, N-1 do begin

				idx = i*N^2 + j*N + k;  

				pos[0, idx] = double(i+0.5) * fac
				pos[1, idx] = double(j+0.5) * fac
				pos[2, idx] = double(k+0.5) * fac

			end
		end
	end

	return, pos
end

; Sample the power spectrum in a sphere in k space, given from the 
; resolution N. We generate two Gaussian random
; numbers using the Box-Mueller method. The FFT of the complex 3D Fourier
; modes has to be real, so the 3D cube of the modes has to be Hermitian. Now
; what does that mean in 3D ? In particular, the k=0 plane becomes rather involved.
; We the scale back the power spectrum to initial redshift using the Zeldovich 
; approximation (Hui & Bertschinger 96). This is known to cause transients, 
; and 2LPT is needed to make it better (Crocce+ 2006).
; If we would start from a glass, we'd need to deconvolve the CIC kernel 
; as well. The displacement fields are resolution independent as well.

function displacement_fields, N, boxsize, Da, showPk=showPk, fout=fout
	
	kmin = 2D*!pi/boxsize
	kmax = N/2D * kmin ; Nyquist mode

	displ = make_array(3, N, N, N, val=0D)
	kmag = make_array(N, N, N, val=0D)

	kdata = make_array(3, N, N, N, val=0D, /complex) ; sampled power spectrum in k-space, 

	kvec = make_array(3, val=0D)

	random_numbers = randomn(14041981, 2, 256, 256, 256, /double) ; res. independent
	
	;readcol, "kspace.dat", ax, ii,jj,kk, km, kx,ky,kz, PowSpec, dplus, fac, aa, pp, dd
 	
	;run = 0L
	
	for comp = 0, 2 do begin

		print, comp

		kdata_rl =  make_array(N, N, N, val=0D) ; keep real and imaginary parts separately
		kdata_im =  make_array(N, N, N, val=0D)

		kmag[*,*,*] = 0D

		for i = 0UL, N-1 do  $
		for j = 0UL, N-1 do  $
		for k = 0UL, N/2-1 do begin 
           
			if (i eq 0) and (j eq 0) and (k eq 0) then $ ;  no DC current
				continue

			if (i eq N/2) or (j eq N/2) or (k eq N/2) then $
				continue

			if i lt N/2. then kvec[0] = i * kmin $ ; k vector at i,j.k
						 else kvec[0] = -(N - i) * kmin
		
			if j lt N/2. then kvec[1] = j * kmin $
						 else kvec[1] = -(N - j) * kmin

			if k lt N/2. then kvec[2] = k * kmin $
					     else kvec[2] = -(N - k) * kmin
	
			kmag[i,j,k]  = sqrt(kvec[0]^2 + kvec[1]^2 + kvec[2]^2)
	
			if kmag[i,j,k] gt kmax then $ ; Only do a sphere in k space
				continue    

;if (k eq 0) and (i eq 0) and (j ge N/2) then continue
;if (k eq 0) and (i ge N/2) then continue

			y1 = random_numbers[0, i,j,k] ; not the Box Mueller transform
			y2 = 2*!pi*random_numbers[1, i,j,k]

			;phase = pp[run]

			Pk  =  Power_Spectrum_EH(kmag[i,j,k])

			Delta =kmin^1.5 * sqrt(Pk/2D) / Da ; scale to a, Zeldovich approx

			; Set kdata so we get a real field after inverse FFT

            kdata_rl[i,j,k] = -kvec[comp] / kmag[i,j,k]^2D * Delta * y1
	    	kdata_im[i,j,k] = kvec[comp] / kmag[i,j,k]^2D * Delta * y1
			
			iconj = N - i  ; conjugated indizes of i j k
			if i eq 0 then $ 
				iconj = 0

			jconj = N - j 
			if j eq 0 then $ 
				jconj = 0

			kconj = N - k 
			if k eq 0 then $
				kconj = 0
			
   	       	kdata_rl[iconj,jconj,kconj] = kdata_rl[i,j,k]
			kdata_im[iconj,jconj,kconj] = -kdata_im[i,j,k]

		;	run++
		end ; for i j k

		kdata[comp,*,*,*] = Complex(kdata_rl, kdata_im, /double)

		data = reform( FFT(kdata[comp, *,*,*], /inverse, /double ) ) 

		bad = where( abs(imaginary(data)) gt 1e-5 * abs(real_part(data)), cnt)

		if cnt ne 0 then $ ; check if we got the symmetries correct
			stop 
	
		displ[comp, *,*,*] = real_part(data)
	end ; comp 

	if keyword_set(showPk) then begin
			
			PK_data = sqrt(kdata[0,*,*,*]^2 + kdata[1,*,*,*]^2 + kdata[2,*,*,*]^2)

			good = where(PK_data ne 0 )

			plot, kmag, PK_data, /xlog, /ylog, xrange=[kmin, kmax], psym=3, $
				xtitle='k [Mpc!U-1!N]', ytitle='P(k, z!Dinit!N)', $
				yrange=minmax(PK_data[good])

			bin_pk = bin_arr(PK_data, pos=kmag, bin_pos=bin_pos, nbins=60, /log)
			oplot, bin_pos, bin_pk, col=color(1)

			dk = alog10(kmax/kmin)/99D
			k_arr = kmin * 10D^(indgen(100)*dk)
			pk_scaled = Power_Spectrum_EH(k_arr)* kmin^1.5/Da
			oplot, k_arr, pk_scaled/(2*!pi^2D), col=color(0)
	end

	return, displ
end

pro set_Pk_vars, Omega_M, Omega_B, hbpar, sigma8, prim_idx, unit_length

	common Pk_vars, pk0, s, agam, Omh, idx

	Mpc2cm = 3.0856802d24 

	idx = prim_idx

	if Omega_B eq 0 then $ ; Pk formulas work only for Omega_B > 0
		stop 

	Omh = Omega_M * hbpar
	Omh2 = Omega_M * hbpar^2
	Obh2 = Omega_B * hbpar^2

	s = 44.5 * alog(9.83D / Omh2) / sqrt(1D + 10D * obh2^0.75) * hbpar ; EH98 eq. 26

	agam = 1D - 0.328 * alog(431D * Omh2) * Omega_B/Omega_M $ ; EH98 eq. 31
			  + 0.38D * alog(22.3D* Omh2) * (Omega_B/Omega_M)^2

	R8 = 8D * Mpc2cm / unit_length ; 8 Mpc/h 
	pk0 = 1D ; normalise with = 1
	pk0 = sigma8^2D / Sigma2_TopHat(R8);
	
	print, "  P(k) Parameters: "
	print, "   prim. idx     = "+strn(idx)
	print, "   s             = "+strn(s) 
	print, "   alpha_Gamma   = "+strn(agam)
	print, "   R8            = "+strn(R8)+" code units"
	print, "   Sigma2_TH(R8) = "+strn(sigma8^2/pk0)
	print, "   Pk0           = "+strn(pk0) 

	return
end

; Eisenstein & Hu 1998, Apj 496; Transfer function for DM and Baryons, no
; neutrinos. 
function Transfer_Function, k 

	common parameters, Tandav,  hbpar,   Omega_M, Omega_L, Tcmb
	common Pk_vars, pk0, s, agam, Omh, idx

	Mpc2cm = 3.0856802d24 

	k_mpc = k * Mpc2cm / Tandav.Length  

	theta = Tcmb/2.7D
		
	gam = Omh*(agam + (1D - agam) / (1D + (0.43D*s*k_mpc)^4D))  ; EH98 eq. 30
	q = k_mpc * theta^2D /gam ; EH98 eq. 28

	L0 = alog(2D*exp(1D) + 1.8D * q) ; EH98 eq. 29
	C0 = 14.2D + 731D / (1D + 62.5D * q) 
	T0 = L0 / (L0 + C0 * q^2D)

	return, T0
end

; Eisenstein & Hu 1998/99 power spectrum at z=0 with n=1. 
; Note that this is P^vel_cb = P(k)/k^2, because we need
; the displacement/velocity. This cancels out with 1D -> 3D
; P^3D(k) = 4 * pi * k^2 * P1D(k) 
function Power_Spectrum_EH, k
	
	common Pk_vars, pk0, s, agam, Omh, idx
	
	T0 = Transfer_Function(k)

	return, pk0 * k^(double(idx)) * T0^2D ; EH99 eq. 25 + 29
end

; EH99 eq.9
function g, a
	
	common parameters, tandav,  hbpar, Omega_M, Omega_L, Tcmb
	
	g2 = Omega_M / a^3 + (1D - Omega_M - Omega_L)/a^2 + Omega_L

	return, sqrt(g2) 
end

; EH99 eq. 10
function growth_integral, a

	common parameters, tandav,  hbpar, Omega_M, Omega_L, Tcmb

	return, ( a / (Omega_M + (1D - Omega_M - Omega_L)*a + Omega_L* a^3))^1.5D
end

function Growth_Function, a
	
	common parameters, tandav, hbpar, Omega_M, Omega_L, Tcmb

	z_eq = 2.5d4 * Omega_M * hbpar^2 * (Tcmb/2.7d)^(-4D) ; EH99 eq.1

	D1_a = Qromb('growth_integral', 0D, a, /double, K=10, Jmax=30D) ; EH99 eq.8

	D1_a *= 5.0/2.0*Omega_M * (1+z_eq) * g(a)

	return, D1_a
end


; find square of mass variance inside a tophat
; window of scale R_Tophat
function sigma_TH_integrant, k 

	common sigma_TH_param, R_Tophat

	kr = R_Tophat * k

	if kr lt 1d-8 then $
		return, 0

	w = 3D * (sin(kr)/kr^3 - cos(kr)/kr^2) ;  EH99, eq. 36
	x = 4D * !pi * k^2 * w^2 * Power_Spectrum_EH(k) ; EH99, eq. 34

	return, x
end

function Sigma2_TopHat, R

	common sigma_TH_param, R_Tophat

	R_Tophat = R ; in Mpc

	kmin = 0D
	kmax = 500D/R

	sigma2 = Qromb('sigma_TH_integrant', kmin, kmax, /double, K=10)

	return, sigma2
end

; variance integrant without window function

function sigma_integrant, k

	return, 4*!pi*k^2 * Power_Spectrum_EH(k) 
end

; find initial redshift from the constraint of linearity
; (no shell crossing) and minimal numerical noise. 
; Limit mass variance in the box to 0.1 at the
; initial scale factor

function find_initial_time, N, boxsize

	common sigma_TH_param, R_Tophat

	max_variance = 0.1D ; in box 

	kmin = 2D*!pi/boxsize
	kmax = N/2D * kmin 

	left = 1e-3
	right = 1D 

	error = 1

	while abs(error) gt 1e-3 do begin ; bisection
	
		a = 0.5 * (right+left)

		Dplus = Growth_Function(1) / Growth_Function(a) 

		int = Qromb('sigma_integrant', kmin, kmax, /double, K=10)
		sigma = sqrt(int / Dplus / (2*!pi^2))

		error = (sigma - max_variance)/max_variance

		if error lt 0 then $
			left = a $
		else $
			right = a
	end 

	return, a
end

; Find the smallest linear mode at z=0. This is
; the minimal boxsize, even though you probably
; want something larger.

function find_min_boxsize, N

	max_variance = 1D ; linearity

	a = 1D ; z = 0

	left = 1D
	right = 1d50 ; code units !

	error = 1

	while abs(error) gt 1e-5 do begin ; bisection
	
		boxsize = 0.5D * (right+left)

		kmin = 2D*!pi/boxsize
		kmax = N/2D * kmin 

		int = Qromb('sigma_integrant', kmin, kmax, /double, K=10)
		sigma = sqrt(int/(2*!pi^2))

		error = (sigma - max_variance)/max_variance
	
		if error gt 0 then $
			left = boxsize $
		else $
			right = boxsize
	end 

	return, boxsize
end

; Cloud in Cell the IDL way, all loops implicit
function CIC, N, pos, ingrid

	grid = reform(ingrid, N, N, N)

	u = reform(pos[0,*], N^3)
    v = reform(pos[1,*], N^3)
    w = reform(pos[2,*], N^3)

    i = long(u)
    j = long(v)
    k = long(w)

	bad = where(i eq N, cnt)
	if cnt gt 0 then $
		i[bad] = N-1

	bad = where(j eq N, cnt)
	if cnt gt 0 then $
		j[bad] = N-1

	bad = where(k eq N, cnt)
	if cnt gt 0 then $
		k[bad] = N-1

	u -= i
	v -= j
	w -= k

	ii = i + 1
	jj = j + 1
	kk = k + 1

	bad = where(ii ge N, cnt)
	if cnt gt 0 then $
		ii[bad] -= N

	bad = where(jj ge N, cnt)
	if cnt gt 0 then $
		jj[bad] -= N

	bad = where(kk ge N, cnt)
	if cnt gt 0 then $
		kk[bad] -= N

    f1 = (1 - u) * (1 - v) * (1 - w);
    f2 = (1 - u) * (1 - v) * w;
    f3 = (1 - u) * v * (1 - w);
    f4 = (1 - u) * v * w;
    f5 = u * (1 - v) * (1 - w);
    f6 = u * (1 - v) * w;
    f7 = u * v * (1 - w);
    f8 = u * v * w;

	result =  f1 * grid[ i, j,k] + f2 * grid[ i, j,kk] $
			+ f3 * grid[ i,jj,k] + f4 * grid[ i,jj,kk] $
			+ f5 * grid[ii, j,k] + f6 * grid[ii, j,kk] $
			+ f7 * grid[ii,jj,k] + f8 * grid[ii,jj,kk];

    return, result
end

pro constrain_particles_to_box, pos, boxsize

	repeat begin 
		
		bad = where(pos lt 0D, cnt1)

		if cnt1 ne 0 then $ 
			pos[bad] += Boxsize

		bad = where(pos gt boxsize, cnt2)

		if cnt2 ne 0 then $ 
			pos[bad] -= Boxsize

		if (cnt1+cnt2 gt 0) then $
			print, "Npart out of box "+strn(cnt1+cnt2)

	end until (cnt1 eq 0 AND cnt2 eq 0)

	return
end

pro test_distribution
	
	common globals, gadget, tandav, cosmo

	fIDL = 'IC_Cosmo_Box_IDL_64'
	print, fIDL, " black"
	velIDL = tandav.readsnap(fidl, "VEL", head=headIDL)

	print, headIDL.npart[1]
	print, max(velIDL)
	print, mean(velIDL[0,*]), mean(velIDL[1,*]), mean(velIDL[2,*])
	print, stddev(velIDL[0,*]), stddev(velIDL[1,*]), stddev(velIDL[2,*])

	fNG = 'IC_Cosmo_Box'
	print, fNG, " green"
	velNG = tandav.readsnap(fng, "VEL", head=headNG)

	print, headNG.npart[1]
	print, max(velNG)
	print, mean(velNG[0,*]), mean(velNG[1,*]), mean(velNG[2,*])
	print, stddev(velNG[0,*]), stddev(velNG[1,*]), stddev(velNG[2,*])

	one = velIDL * 0+1D

	tmp = bin_arr(one, pos=velIDL[0,*], bin_pos=bin_pos, nbins=200, cnt=cntIDL0)
	plot, bin_pos, double(cntIDL0)/headIDL.npart[1], psym=10, yrange=[0, 0.03d], $
		xrange=[-2000D,2000]

	tmp = bin_arr(one, pos=velIDL[1,*], bin_pos=bin_pos, nbins=200, cnt=cntIDL1)
	oplot, bin_pos, double(cntIDL1)/headIDL.npart[1], psym=10

	tmp = bin_arr(one, pos=velIDL[2,*], bin_pos=bin_pos, nbins=200, cnt=cntIDL2)
	oplot, bin_pos, double(cntIDL2)/headIDL.npart[1], psym=10

	tmp = bin_arr(one, pos=velNG[0,*], bin_pos=bin_pos, nbins=200, cnt=cntNG0)
	oplot, bin_pos, double(cntNG0)/headNG.npart[1], psym=10, col=color(0)

	tmp = bin_arr(one, pos=velNG[1,*], bin_pos=bin_pos, nbins=200, cnt=cntNG1)
	oplot, bin_pos, double(cntNG1)/headNG.npart[1], psym=10, col=color(0)

	tmp = bin_arr(one, pos=velNG[2,*], bin_pos=bin_pos, nbins=200, cnt=cntNG2)
	oplot, bin_pos, double(cntNG2)/headNG.npart[1], psym=10, col=color(0)

	;plot, length(velIDL)
	;oplot, length(velng), col=color(0)

	return
end
