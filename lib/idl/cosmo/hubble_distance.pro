function hubble_distance, H0  ; hubble distance
    @set_cgs  
    return, c/double(H0)
end
