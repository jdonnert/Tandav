function scale_factor, t

    print, "Not implemented"
    return, -1
end

function lookback_time, z, self, inv=inv
    return, 1
end

function hubble_param, z, H0
    
    return, H0 / E_peebles(z)   ; 1993, pp. 310 - 321
end

function E_peebles, z, Omega_M, Omega_k, Omega_L
???
???
