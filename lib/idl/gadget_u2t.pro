function gadget_U2T, U, inv=inv, xH=xH, uvel=uvel, gamma=gamma, radiative=radiative

    @set_cgs

    if not exist(U) then begin
        print, "U2T, U, inv=inv, xH=xH, uvel=uvel, gamma=gamma, rad=rad, uVel=uVel"
        return, -1
    end

    if keyword_set(rad) then begin
        print, "Radiative Runs not implemented yet"
		return, -1
    end
    
    if not exist(uvel) then $
        uvel = 1d5

    if not keyword_set(gamma) then $
        gamma = 5.0/3.0

    if not keyword_set(xH) then $
        xH = 0.76

  	yhelium = ( 1. - xH ) / ( 4 * xH )

    mean_mol_weight = (1. + 4. * yhelium) / (1. + 3. * yhelium + 1)

    U2T_factor = (gamma-1) * uvel^2 * mp*mean_mol_weight/k_boltz

    if  keyword_set(inv) then $
        U2T_factor  = 1D / U2T_factor

    return,  U * U2T_factor 
end
