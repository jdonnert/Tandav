; P-Tandav Class
pro TandavCodeObject__define

    void = {TandavCodeObject,       $   
                ; units
                name    : '  ',     $
                length  : double(0),$
                mass    : double(0),$
                velocity: double(0),$
                time    : double(0),$
                energy  : double(0),$   
                ; constants 
                fr      : double(0),$   ; hydrogen mass fraction
                gamma   : double(0),$   ; adiabatic index 
                xH      : double(0),$   ; hydrogen fraction
                inherits IDL_Object $   ; for public variables
    }

    return
end

function TandavCodeObject::INIT
    
    self.Set, 2

    return, 1
end

; expose object variables
pro TandavCodeObject::GetProperty, length=length, mass=mass,$
        velocity=velocity, time=time, energy=energy, fr=fr,$
        gamma=gamma, xH=xH

    ; Units
    if arg_present(length)  then length = self.length
    if arg_present(mass)    then mass = self.mass
    if arg_present(velocity)then velocity = self.velocity
    if arg_present(time)    then time = self.time
    if arg_present(energy)  then energy = self.energy
    if arg_present(fr)      then fr = self.fr
    if arg_present(gamma)   then gamma = self.gamma
    if arg_present(xH)      then xH = self.xH

    return
end

pro TandavCodeObject::SetProperty, length=length, mass=mass,$
        velocity=velocity, time=time, energy=energy, fr=fr,$
        gamma=gamma, xH=xH, preset=preset

    ; Units
    if keyword_set(length)  then self.length = length
    if keyword_set(mass)    then self.mass = mass
    if keyword_set(velocity)then self.velocity = velocity
    if keyword_set(time)    then self.time = time
    if keyword_set(energy)  then self.energy = energy
    if keyword_set(fr)      then self.fr = fr
    if keyword_set(gamma)   then self.gamma = gamma
    if keyword_set(xH)      then self.xH = xH

    self.set, 0, /silent   ; just recompute derived values

    return
end


; show current values
pro TandavCodeObject::Show
    @set_cgs

    print, self.name+' :'
    print, '    length   = '+strn(self.length, format='(1e8.2)')+' cm'
    print, '    mass     = '+strn(self.mass, format='(1e8.2)')+' g'
    print, '    velocity = '+strn(self.velocity, format='(1e8.2)')+' cm/s'
    print, '    xH       = '+strn(self.xH)
    print, '    adiab.idx= '+strn(self.gamma)

    return
end

; set different systems
pro TandavCodeObject::Set, val, silent=silent
    
    @set_cgs

    if n_params() lt 1 then begin
        print, "List of Tandav unit systems"
        print, "    0 - l=kpc/h, m=10^10 Msol, v=1km/sec, xH=0.76, gam=5/3"
    end

    case val of
        0 : begin
            ; only recompute derived quantities
            end
        1 : begin
                self.name       = 'Cgs Units'
                self.length     = 1         
                self.mass       = 1
                self.velocity   = 1

                self.fr         = 0.76D
                self.gamma      = 5D/3D             ; monoatomic gas
                self.xH         = 0.76D
            end
        2 : begin
                self.name       = 'Cosmo Units, kpc, 1e10 Msol, km/s'
                self.length     = kpc2cm            ; 1.0 kpc /h
                self.mass       = 1d10*Msol         ; 10^10 solar masses
                self.velocity   = 1d5               ; 1 km/sec

                self.fr         = 0.76D
                self.gamma      = 5D/3D             ; monoatomic gas
                self.xH         = 0.76D
            end
        
    endcase
    
    self.time   = self.length/self.velocity
    self.energy = self.mass*self.length^2 / self.time^2

    if not keyword_set(silent) then $
        self.show
    
    return
end

;compute physical density
function TandavCodeObject::Density, rho, h=h, z=z, electrons=electrons
    
   return, tandav_density(rho, h=h, z=z, electrons=electrons, xH=self.xH, $
       uMass=self.mass, uLength=self.length)
end

; Simple Thermodynamics
function TandavCodeObject::U2T, U, inv=inv, radiative=radiative

    return, tandav_U2T( U, inv=inv, xH=self.xH, uvel=self.velocity, $
        gamma=self.gamma, radiative=radiative )
end

function TandavCodeObject::T2U, T, inv=inv

    return, self.U2T( T, inv=1 )
end 

function TandavCodeObject::ThermalEnergyDensity, rho, T, h=h, z=z, TisU=TisU

    return, tandav_thermal_energy_density( rho, T, self.xH, h=h, z=z, TisU=TisU )
end

function TandavCodeObject::Pressure, rho, U, z=z, h=h

    return, tandav_pressure(rho, u, gamma=self.gamma, xH=self.xH, z=z, h=h,$
        uVel=self.Velocity, uMass=self.Mass, uLength=self.Length)
end

; Include snapshot reading
function TandavCodeObject::ReadSnap, fname, block, head=head, $
    parttype=parttype, debug=debug

    return, read_tandav_snapshot(fname, block, head=head, $
        						parttype=parttype,debug=debug)
end

; Include snapshot writing
function TandavCodeObject::MakeHead

    return, make_head()
end

pro TandavCodeObject::WriteHead, fname, head, swap_endian=swap_endian

    write_head, fname, head, swap_endian=swap_endian 

    return
end

pro TandavCodeObject::AddBlock, fname, data, blockid, debug=debug, $
    swap_endian=swap_endian

    add_block, fname, data, blockid, debug=debug, swap_endian=swap_endian

    return
end
