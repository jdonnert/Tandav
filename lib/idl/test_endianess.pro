function test_endianess, fname, debug=debug

    four_bytes = long(0)    

    openr, fp, fname, /get_lun
    readu, fp, four_bytes
    close, fp
    free_lun, fp, /force
        
    if four_bytes eq 8 then begin 
        if keyword_set(debug) then $
            print, fname+' is same endianess'

        return, 0   
    end
        
    if four_bytes eq 134217728 then begin 
        if keyword_set(debug) then $
            print, fname+' is other endianess'

        return, 1
    end

    print, 'Error, Snapshot corrupt or format 1'

    return, -1
end
