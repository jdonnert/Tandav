; make simple DM only cosmological initial conditions
; Eisenstein & Hu power spectrum, similar to N-GenIC

pro make_ICs

	common parameters hbpar, Omega_M, Omega_Baryon

	tandav = obj_new('TANDAVCODEOBJECT')

	tandav.set, 1 ; cosmo units

	N = 128L
	boxsize = 150000D

	z = 63D
	time = 1/ (1 + z)

	G = tandav.Grav ; tandav units
	sigma8 = 0.8D
	hbpar = 0.7D
	hubble = 3.2407766d-18 * hbpar * tandav.time
	Omega_M = 0.3D
	Omega_L = 0.7D
	Omega_B = 0D
	Omega_0 = Omega_M + Omega_L + Omega_B
		
	npart = N^3
	cellsize = boxsize/N

	pos = make_mesh(N)

	displ = displacement_fields()

	displ_vec = make_array(3, npart, val=0D)
	displ_vec[0,*] = CIC(N, pos/cellsize, displ[0,*,*,*])
	displ_vec[1,*] = CIC(N, pos/cellsize, displ[1,*,*,*])
	displ_vec[2,*] = CIC(N, pos/cellsize, displ[2,*,*,*])

	vel = make_array(3, npart, val=0D)

	for i = 0, npart-1 do begin

		pos[0,i] += binned_displs[0,i]
		pos[1,i] += binned_displs[1,i]
		pos[2,i] += binned_displs[2,i]

		vel[0,i] = binned_displs[0,i] * displ2vel
		vel[1,i] = binned_displs[1,i] * displ2vel
		vel[2,i] = binned_displs[2,i] * displ2vel

	end

	id = ulindgen(npart)+1

	head = tandav.make_head()
	head.npart = [0,1,0,0,0,0] * npart
	head.massarr = [0,1,0,0,0,0] * total_mass/npart
	head.time = time
	head.redshift = z
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

				pos[0, idx] = double(i) * fac
				pos[1, idx] = double(j) * fac
				pos[2, idx] = double(k) * fac

			end
		end
	end

	return, pos
end

; Construct the displacement fields. Sample the power spectrum in a sphere
; in k space, given from the resolution N. The FFT of the complex 3D Fourier
; modes has to be real, so the 3D cube of the modes has to be Hermitian. Now
; what does that mean in 3D ? In particular, the i=0 plane becomes rather involved.
; Scale back the power spectrum to initial redshift using the Zeldovich approximation. 

function displacement_fields, N, boxsize, show_Pk=show_Pk

	seeed = 140419081 ; Dickes B ! Oben an der Spree ! Im Sommer tust du gut, ...

	kmin = 2D*!pi/boxsize
	kmax = !pi*N/boxsize

	displ = make_array(N, N, N, val=0D)
	cdata = make_array(N, N, N, val=0D,/complex) ; sampled power spectrum in k-space
	cdata_rl =  make_array(N, N, N, val=0D)
	cdata_im =  make_array(N, N, N, val=0D)

	kvec = make_array(3, val=0D)

	for comp = 0, 2 do begin

		for i = 0, N/2 do $
		for j = 0, N-1 do $
		for k = 0, N-1 do begin ; Generate k values first
            
			if i eq N/2 or j eq N/2 or k eq N/2 then $
				continue

			if i eq 0 and j eq 0 and k eq 0 then $ ;  no DC current
				continue

			if i ne 0 then iconj = N - i $ ; Define conjugated indizes of i j k
					  else iconj = 0
			if j ne 0 then jconj = N - j $
					  else jconj = 0
			if k ne 0 then kconj = N - k $
					  else kconj = 0
            
			if i LE N/2. then kvec[0] = i * kmin $ ; Define grid
						 else kvec[0] = -iconj * kmin
	
			if j LE N/2. then kvec[1] = j * kmin $
						 else kvec[1] = -jconj * kmin
	
			if k LE N/2. then kvec[2] = k * kmin $
					     else kvec[2] = -kconj * kmin

			kmag[i,j,k]  = sqrt(kvec[0]^2 + kvec[1]^2 + kvec[2]^2)
	
			if kmag[i,j,k] GT kmax then $ ; Only do a sphere in k space
				continue    

			A = randomu(seeed)
			phase = randomu(seeed)

			Pk  = -alog(A) * Power_Spectrum_EH(k)
	
			dPk = nrm * sqrt(Pk) / Da ; scales the spectrum back to z_init

			; Set power so we get a real field after inverse FFT

            if i gt 0 then begin    ; grid is hermitian in i > ngrid/2

                cdata_rl[i,j,k] = -kvec[comp] / kmag[i,j,k]^2 * dPk * sin(phase)
		    	cdata_im[i,j,k] = kvec[comp] / kmag[i,j,k]^2 * dPk * cos(phase)

                cdata_rl[iconj,jconj,kconj] = cdata_rl[i,j,k] ; DOUBTFUL, PLOT !!
    			cdata_im[iconj,jconj,kconj] = -1*cdata_im[i,j,k]

            end else begin  ; i = 0 plane is where we have to get the symmetry right

                if j eq 0 then begin ; first row
                    
                    if k gt N/2. then $
                        continue
                
                    cdata_rl[i,j,k] = -kvec[comp] / kmag[i,j,k]^2 * dPk * sin(phase)
		    	    cdata_im[i,j,k] = kvec[comp] / kmag[i,j,k]^2 * dPk * cos(phase)
                
                    cdata_rl[i,j,kconj] = cdata_rl[i,j,k]
    			    cdata_im[i,j,kconj] = -1*cdata_im[i,j,k]

                end else begin  ; j != 0 here

                    if j gt N/2. then $     ; rest of the plane
                        continue
                
                    cdata_rl[i,j,k] = -kvec[comp] / kmag[i,j,k]^2 * dPk * sin(phase)
		    	    cdata_im[i,j,k] = kvec[comp] / kmag[i,j,k]^2 * dPk * sin(phase)

                    cdata_rl[i,jconj,kconj] = cdata_rl[i,j,k]
    		    	cdata_im[i,jconj,kconj] = -1*cdata_im[i,j,k]
                end
            end ; if i
		end ; ijk

		cdata[comp,*,*,*] = Complex(cdata_rl, cdata_im)
		data = reform( FFT(cdata[0,*,*,*], /inverse, /double ) ) 

		
        for i = 0, ngrid^3-1 do $ ; check if we got the symmetries correct
            if abs(imaginary(data[i])) gt 1e-7 * abs(real_part(data[i])) then $ 
                stop
	
		displ[comp, *,*,*] = real_part(data) / N^3

		if keyword_set(show_Pk) then begin



		end

	end ; comp 

	return, displ
end

function power_spectrum_EH, k
	
	common parameters hbpar, Omega_M, Omega_Baryon

	pk0 = 1
	tk = 1

	Pk = pk0 * k * tk^2

	return, Pk
end


; Cloud in Cell the IDL way
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

	bad = where(i ge N, cnt)

	if cnt gt 0 then $
		i[bad] -= N

	bad = where(j ge N, cnt)
	
	if cnt gt 0 then $
		j[bad] -= N

	bad = where(k ege N, cnt)
	
	if cnt gt 0 then $
		k[bad] -= N

    f1 = (1 - u) * (1 - v) * (1 - w);
    f2 = (1 - u) * (1 - v) * w;
    f3 = (1 - u) * v * (1 - w);
    f4 = (1 - u) * v * w;
    f5 = u * (1 - v) * (1 - w);
    f6 = u * (1 - v) * w;
    f7 = u * v * (1 - w);
    f8 = u * v * w;

	result = 1;

    return, grid[i,j,k]
end
