pro add_block,fname,data,blockid,debug=debug,SWAP_ENDIAN=SWAP_ENDIAN, force=force

    type = size(data,/type)
    if type ne 4 and type ne 13 then $
        print, 'WARNING : Input data is neither float nor ulong !'
    
    siz_data = size(data)

    if siz_data[0] eq 2 then $  ; vector
        nBytes = siz_data[1]*siz_data[2]*4L + 8L $
    else $                      ; scalar
        nBytes = siz_data[1] * 4L + 8L

    b4n = byte(blockid+"    ")
    b4n = b4n(0:3)

    if keyword_set(debug) then $
        print, nBytes, ':', string(b4n)

    pos = long64(0)
    b4 = bytarr(4)
    bl = long(0)

    swap_endian = test_endianess(fname, debug=debug)

    openu, fd, fname, /f77, /get_lun, swap_endian=swap_endian

    if keyword_set(debug) then begin
        while not (eof(fd)) do begin
            readu, fd, b4, bl
         
            thislabel = string(b4)
          
            print,thislabel," <=> ", string(b4n), bl

            point_lun,-fd,pos

            point_lun,fd,long64(pos)+long64(bl)
        end
    end else begin
        skip_lun, fd, /EOF
    end

    writeu, fd, b4n, nBytes
    writeu, fd, data
    close, fd

   free_lun,fd

   return
END   
