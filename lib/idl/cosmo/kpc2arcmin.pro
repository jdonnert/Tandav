function  arcmin2kpc, X, z, H0, Omega_M, Omega_L, inv=inv
    
    if n_params() lt 5 $
    or ~keyword_set(X) or ~keyword_set(z) or ~keyword_set(H0) $
    or ~keyword_set(Omega_M) or ~keyword_set(Omega_L) then begin
        print, 'N_kpc = arcmin2kpc(N_arcmin, z, H0, Omega_M, Omega_L, inv=inv)'
        return, -1
    end

    @set_cgs

    d_ang = angular_diameter_distance(z, H0, Omega_M, Omega_L)

    arcmin2rad = 2.*!pi/360./60.

    if not keyword_set(inv) then $  ; X = N_arcmin
        return, X * d_ang/kpc2cm * arcmin2rad $
    else $                          ; X = N_kpc
        return, X / (d_ang/kpc2cm * arcmin2rad)
end

function kpc2arcmin, X, z, H0, Omega_M, Omega_L
    return, arcmin2kpc(X, z, H0, Omega_M, Omega_L, inv=1)
end
