; convert density in gadget units to physical
function gadget_density, rho, h=h, z=z, electrons=electrons, xH=xH, uMass=uMass, uLength=uLength

    @set_cgs

    if not (exist(rho) or exist(xH) or exist(uMass) or exist(uLength) ) then begin
        print, "gadget_density( rho, h=h, z=z, electrons=electrons, xH=xH, "
        print, "                 uMass=uMass, uLength=uLength )"
        return, -1
    end
        
    if not keyword_set(z) then $
        z = 0

    if not keyword_set(h) then $
        h = 0.70    ;fall back

    ; conversion n_pat -> n_electrons
    n2ne = (xH+0.5*(1-xH))/(2*xH+0.75*(1-xH)) 

    ; mean molecular weight in hydr. mass
    umu = 4./(5.*xH+3.) 

    fac = 1
    if keyword_set(electrons)  then $
        fac = n2ne/(umu*mp)

    return, double(rho) * (1+z)^3*uMass/uLength^3 *h^2 * fac
end
