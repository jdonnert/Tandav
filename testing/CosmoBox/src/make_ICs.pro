; make simple DM only cosmological initial conditions
; Eisenstein & Hu power spectrum, Zeldovich Approximation.
; [gadget]: use gadget units in velocity and hbpar=1
; [showPk]: plot P(k) from sampled data and analytic formula

pro make_ICs, gadget=gadget, showPk=showPk

	common parameters, tandav,  hbpar,  Omega_M, Omega_L, Tcmb

	tandav = obj_new('TANDAVCODEOBJECT')
	
	Grav = tandav.Grav ; tandav units

	Mpc2cm = 3.0856802d24 

	N = 64L
	boxsize = 150000D
	npart = N^3

	print, "npart = "+strn(npart)

	z_init = 63D
	a_init = 1D/(1D + z_init)

	print, "initial a = "+strn(a_init)

	prim_idx = 1D  ; index of primordial power spectrum
	sigma8 = 0.8D  ; normalisation of P(k) at 8 Mpc
	Tcmb = 2.728D  ; temperature of the CMB [K]
	hbpar = 0.7D   ; Hubble parameter
	Omega_L = 0.7D ; Dark Energy at z=0
	Omega_M = 0.3D ; Total Matter, EH98 call this Omega_0
	Omega_B = 0.04 ; Baryons 
	
	H0 = 100D * 1d5 / Mpc2cm  * tandav.time ; tandav Units

	if not keyword_set(gadget) then $ ; hbpar != 1
		H0 *= hbpar

	rho_crit = 3D * H0^2/(8D*!pi*Grav) ; comoving tandav units
	total_mass = Omega_M * rho_crit * boxsize^3

	print, "H0 = "+strn(H0)
	print, "rho_crit(a_init) = "+strn(rho_crit)
	print, "Mtot = "+strn(total_mass)
	print, "Mpart = "+strn(total_mass/npart)

	set_Pk_vars, Omega_M, Omega_B, hbpar, sigma8, prim_idx, tandav.length

	Da = Growth_Function(1D) / Growth_Function(a_init)

	print, "Da = "+strn(Da) ; growth factor for ZA

	pos = make_mesh(N, Boxsize)

	displ_grid = displacement_fields(N, Boxsize, Da, showPk=showPk)

	cellsize = boxsize/N

	displ = make_array(3, npart, val=0D)
	displ[0,*] = CIC(N, pos/cellsize, displ_grid[0,*,*,*])
	displ[1,*] = CIC(N, pos/cellsize, displ_grid[1,*,*,*])
	displ[2,*] = CIC(N, pos/cellsize, displ_grid[2,*,*,*])

	pos += displ

	print, "mean displ = "+strn(mean(abs(displ_grid)))+" kpc/h"
	print, "max(displ) = "+strn(max(displ))+" kpc/h" 

	f0 = ( Omega_M / a_init^3 / g(a_init)^2)^(0.6D)  ; EH98, eq. 28
	displ2vel = a_init * g(a_init) * H0 * f0; EH99, eq. 29
	
	if keyword_set(gadget) then $
		vel2comov = 1D/sqrt(a_init) $ ; Gadget comoving vel ~ 1/sqrt(a)
	else $
		vel2comov = 1D/a_init^2 ; tandav comoving vel ~ 1/a^2

	vel = displ * displ2vel * vel2comov

	print, "F0(a_init) = "+strn(f0) 
	print, "displ2vel = "+strn(displ2vel)
	print, "vel2comov = "+strn(vel2comov)

	constrain_particles_to_box, pos, boxsize

	id = ulindgen(npart)+1

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

	fout = "IC_Cosmo_Box_IDL"

	print, fout

	tandav.write_head, fout, head
	tandav.add_block, fout, float(pos), "POS"
	tandav.add_block, fout, float(vel), "VEL"
	tandav.add_block, fout, ulong(id), "ID"

	return
end

function make_mesh, N,  boxsize
	
	pos = make_array(3, N^3, val=0D)

	fac = boxsize/N

	for i = 0, N-1 do begin
		for j = 0, N-1 do begin
			for k = 0, N-1 do begin

				idx = i*N^2 + j*N + k;  

				pos[0, idx] = double(i+0.5) * fac
				pos[1, idx] = double(j+0.5) * fac
				pos[2, idx] = double(k+0.5) * fac

			end
		end
	end

	return, pos
end

; Construct the displacement fields. Sample the power spectrum in a sphere
; in k space, given from the resolution N. We generate two Gaussian random
; number using the Box-Mueller method. The FFT of the complex 3D Fourier
; modes has to be real, so the 3D cube of the modes has to be Hermitian. Now
; what does that mean in 3D ? In particular, the i=0 plane becomes rather involved.
; To complicate things, this depends on the FFT used, so the symmetries in IDL
; are not the same as in FFTW, for example. 
; We the scale back the power spectrum to initial redshift using the Zeldovich 
; approximation (Hui & Bertschinger 96). This is known to cause transients, 
; and 2LPT is needed to make it better (Crocce+ 2006).

function displacement_fields, N, boxsize, Da, showPk=showPk
	
	kmin = 2D*!pi/boxsize
	kmax = N/2D * kmin ; Nyquist mode

	displ = make_array(3, N, N, N, val=0D)
	kmag = make_array(N, N, N, val=0D)

	cdata = make_array(3, N, N, N, val=0D, /complex) ; sampled power spectrum in k-space, 
	cdata_rl =  make_array(N, N, N, val=0D) ; keep real and imaginary parts seperately
	cdata_im =  make_array(N, N, N, val=0D)

	kvec = make_array(3, val=0D)
	
	random_field = randomu(14041981L, 2, 3, N, N, N, /double)

	for comp = 0, 2 do begin

		print, comp

		for i = 0, N/2 do $ ; FORTRAN convention
		for j = 0, N-1 do $
		for k = 0, N-1 do begin 
            
			if (i eq N/2) or (j eq N/2) or (k eq N/2) then $
				continue

			if (i eq 0) and (j eq 0) and (k eq 0) then $ ;  no DC current
				continue

			if i ne 0 then iconj = N - i $ ; conjugated indizes of i j k
					  else iconj = 0
			if j ne 0 then jconj = N - j $
					  else jconj = 0
			if k ne 0 then kconj = N - k $
					  else kconj = 0
            
			if i lt N/2. then kvec[0] = i * kmin $ ; k vector at i,j.k
						 else kvec[0] = -iconj * kmin
	
			if j lt N/2. then kvec[1] = j * kmin $
						 else kvec[1] = -jconj * kmin
	
			if k lt N/2. then kvec[2] = k * kmin $
					     else kvec[2] = -kconj * kmin

			kmag[i,j,k]  = sqrt(kvec[0]^2 + kvec[1]^2 + kvec[2]^2)

			if kmag[i,j,k] gt kmax then $ ; Only do a sphere in k space
				continue    

			A = random_field[0, comp, i,j,k]
			phase = 2*!pi * random_field[1, comp, i,j,k]

			Pk  = -alog(A) * Power_Spectrum_EH(kmag[i,j,k]) ; Box Mueller method
	
			Pk_ainit = kmin^1.5 * sqrt(Pk) / Da ; ZA, Mo+ 2010, eq. 4.101

			; Set cdata so we get a real field after inverse FFT

            if i gt 0 then begin ; grid is hermitian in i > ngrid/2

                cdata_rl[i,j,k] = -kvec[comp] / kmag[i,j,k]^2 * Pk_ainit * sin(phase)
		    	cdata_im[i,j,k] = kvec[comp] / kmag[i,j,k]^2 * Pk_ainit * cos(phase)

                cdata_rl[iconj,jconj,kconj] = cdata_rl[i,j,k] ; this is not stored in FFTW
		    	cdata_im[iconj,jconj,kconj] = -1*cdata_im[i,j,k]

            end else begin  ; i = 0 plane is where we have to get the symmetry right

                if j eq 0 then begin ; first row
                    
                    if k gt N/2. then $ 
                         continue ; set already via conjugated indices
                
                    cdata_rl[i,j,k] = -kvec[comp] / kmag[i,j,k]^2 * Pk_ainit * sin(phase)
		    	    cdata_im[i,j,k] = kvec[comp] / kmag[i,j,k]^2 * Pk_ainit * cos(phase)
                
                    cdata_rl[i,j,kconj] = cdata_rl[i,j,k]
    			    cdata_im[i,j,kconj] = -1*cdata_im[i,j,k]

                end else begin  ; j != 0

                    if j ge N/2. then $ 
                        continue ; set already via conjugated indices
                
                    cdata_rl[i,j,k] = -kvec[comp] / kmag[i,j,k]^2 * Pk_ainit * sin(phase)
		    	    cdata_im[i,j,k] = kvec[comp] / kmag[i,j,k]^2 * Pk_ainit * sin(phase)

                    cdata_rl[i,jconj,kconj] = cdata_rl[i,j,k]
    		    	cdata_im[i,jconj,kconj] = -1*cdata_im[i,j,k]
                end
            end ; if i
		end ; for i j k

		cdata[comp,*,*,*] = Complex(cdata_rl, cdata_im)
		data = reform( FFT(cdata[comp,*,*,*], /inverse, /double ) ) 
		
		bad = where( abs(imaginary(data)) gt 1e-5 * abs(real_part(data)), cnt)
		if cnt ne 0 then $ ; check if we got the symmetries correct
			stop 
	
		displ[comp, *,*,*] = real_part(data)

	end ; comp 
	
	if keyword_set(showPk) then begin
			
			PK_data = sqrt(cdata[0,*,*,*]^2 + cdata[1,*,*,*]^2 + cdata[2,*,*,*]^2)

			good = where(PK_data ne 0 )

			plot, kmag, PK_data, /xlog, /ylog, xrange=[kmin, kmax], psym=3, $
				xtitle='k [Mpc!U-1!N]', ytitle='P(k, z!Dinit!N)', yrange=minmax(PK_data[good])

			bin_pk = bin_arr(PK_data, pos=kmag, bin_pos=bin_pos, nbins=60, /log)
			oplot, bin_pos, bin_pk, col=color(1)

			oplot, kmag, Power_Spectrum_EH(kmag)* kmin^1.5/Da, col=color(0)
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
			  + 0.38D * alog(22.3D*Omh2) * (Omega_B/Omega_M)^2

	R8 = 8D * Mpc2cm / unit_length ; 8 Mpc/h 
	pk0 = 1D ; normalise with = 1
	pk0 = sigma8^2D / Sigma2_TopHat(R8);
	
	print, "s = "+strn(s) 
	print, "alpha_Gamma = "+strn(agam)
	print, "R8 = "+strn(R8)+" Mpc"
	print, "Sigma2_TH(R8) = "+strn(sigma8^2/pk0)
	print, "Pk0 = "+strn(pk0) 

	return
end

; Eisenstein & Hu 1998, Apj 496; Transfer function for DM and Baryons, no
; neutrinos. 
function Transfer_Function, k 

	common parameters, tandav,  hbpar,   Omega_M, Omega_L, Tcmb
	common Pk_vars, pk0, s, agam, Omh, idx

	Mpc2cm = 3.0856802d24 

	k_mpc = k * Mpc2cm / Tandav.Length  

	theta = Tcmb/2.7D
		
	gam = Omh*(agam + (1D - agam) / (1D + (0.43D*s*k_mpc)^4D))  ; EH98 eq. 30
	q = k_mpc * theta^2D /gam ; EH98 eq. 28

	L0 = alog(2D*exp(1D) + 1.8D * q) ; EH98 eq. 29
	C0 = 14.2D + 731D / (1D + 62.5D * q) 
	T0 = L0 / (L0 + C0 * q^2)

	return, T0
end

; Eisenstein & Hu 1998/99 power spectrum at z=0 with n=1. 
; Note that this is P^vel_cb = P(k)/k^2, because we need
; the displacement/velocity. This cancels out with 1D -> 3D
; P^3D(k) = 4 * pi * k^2 * P^1D(k) 
function Power_Spectrum_EH, k
	
	common Pk_vars, pk0, s, agam, Omh, idx
	
	T0 = Transfer_Function(k)

	return, pk0 * k^idx * T0^2 ; EH99 eq. 25 + 29
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

	D1_a = Qromb('growth_integral', 0, a, /double, K=10, Jmax=30D) ; EH99 eq.8

	D1_a *= 5.0/2.0*Omega_M * (1+z_eq) * g(a)

	return, D1_a
end

 
function sigma_integral, k

	common sigma_param, R_Tophat

	kr = R_Tophat * k

	if kr lt 1d-8 then $
		return, 0

	w = 3D * (sin(kr)/kr^3 - cos(kr)/kr^2) ;  EH99, eq. 36
	x = 4D * !pi * k^2 * w^2 * Power_Spectrum_EH(k) ; EH99, eq. 34

	return, x
end

function Sigma2_TopHat, R

	common sigma_param, R_Tophat

	R_Tophat = R

	kmin = 0D
	kmax = 500D/R

	sigma2 = Qromb('sigma_integral', kmin, kmax, /double, K=10)

	return, sigma2
end


; Cloud in Cell the IDL way, all loops implicit
function CIC, N, pos, ingrid

	grid = reform(ingrid, N, N, N)

	npart = N^3

	u = reform(pos[0,*], npart)
    v = reform(pos[1,*], npart)
    w = reform(pos[2,*], npart)

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
