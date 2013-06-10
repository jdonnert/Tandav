; P-Gadget snapshot format 2
function get_nfiles, fname
    
    name = file_search(fname, count=count)
    
    if count eq 1 then $ ; try single file
        return, count
    
    for i=0,90000L do begin  ; multiple files
        name = file_search(fname+'.'+strn(i), count=count)

        if count eq 0 then $
            return, i
    end  
    
    return, 0  ; if you have more than 90001 then go away
end

function get_next_block_label, fp, blocksize=blocksize
    
    label = '    ' 
    blocksize = 0UL

    readu, fp, label, blocksize

    return, label
end 

function read_gadget_header, fp

    header = make_head()

    label = get_next_block_label(fp)
    if label ne "HEAD" then $
        print, "Error: First block is not HEAD but "+label

    readu, fp, header
    
    return, header
end
    
function generate_mass_from_header, blockinfo, header

    npart = [0ULL, total(ulong64(header.npart), /cumulative)]

    mass = make_array(npart[6], val=0)

    for i=0,5 do $
        mass[npart[i]:npart[i+1]] = header.massarr[i]
    
    return, mass
end

; define blocks to allow special treatment
function get_Block_info, blockname, npart, want_type=want_type,debug=debug
    
    npart = ulong64(npart) 

    blockinfo = {   name    : blockname,  $
                    ptype   : ulon64arr(6),  $  
                    rank    : ulong64(0), $
                    npart   : ulong64(0), $ 
                    IDLtype : fix(0),     $ ; type code as size()
                    nBytes  : ulong64(0), $ ; per particle
                    known   : 1           $
                }

    switch blockname of
        ; 3D  float all types
        "POS ": 
        "VEL ": begin
                    blockinfo.ptype[*] = 1
                    blockinfo.rank = 3
                    blockinfo.IDLtype = 4
                    break
                end
        ; 1D float all types
        "TSTP": 
        "MASS": begin
                    blockinfo.ptype[*] = 1
                    blockinfo.rank = 1
                    blockinfo.IDLtype = 4
                    break
                end
        ; 1D uint, all types
        "ID  ": begin 
                    blockinfo.ptype[*] = 1
                    blockinfo.rank = 1
                    blockinfo.IDLtype = 13
                    break
                end
        ; 1D uint, type 0 only
        "TNGB": begin 
                    blockinfo.ptype[0] = 1
                    blockinfo.rank = 1
                    blockinfo.IDLtype = 13   
                    break
                end
        ;3D float, type 0 only
        "VBLK": 
        "SHNR":        
        "BFLD": begin
                    blockinfo.ptype[0] = 1
                    blockinfo.rank = 3
                    blockinfo.IDLtype = 4
                    break
                end 
        ; std block is 1D float, type 0 only
        else  : begin   
                    blockinfo.ptype[0] = 1
                    blockinfo.rank = 1
                    blockinfo.IDLtype = 4
                    blockinfo.known = 0
                    break
                end
    end

    case blockinfo.IDLtype of
        4   : blockinfo.nBytes = 4  ; sizeof(float)
        13  : blockinfo.nBytes = 4  ; sizeof(uint)
    end
    
    blockinfo.npart = total(npart*blockinfo.ptype)

    if blockinfo.known eq 0 and keyword_set(debug) then $
        print, 'Assuming unknown block "'+blockname+'" is 1D, float, SPH only'

    return, blockinfo
end

pro skip_block, fp, blocksize
    
    point_lun,-fp, curr_offset      ;get current offset

    point_lun, fp, ulong64(curr_offset) + ulong64(blocksize)
    
    return
end

function read_block, fp, blockinfo

    data = make_array(blockinfo.rank, blockinfo.npart, type=blockinfo.IDLtype)
    
    readu, fp, data

    return, data
end

; sort new data into all data and keep type order
; here rank means length of vector entries not dimensions
pro sort_particles, new_data, new_npart, all_data=all_data, $
    all_npart=all_npart, want_type=want_type, parttotal=parttotal
  
    type = size(new_data,/type)
    nDim = size(new_data,/n_dimension)

    if nDim eq 1 then begin  ; treat scalar block & one particle case
        rank = 1
        nData = (size(new_data))[1] 
    end else begin
        rank = (size(new_data))[1]
        nData = (size(new_data))[2]
    end

    good = where(want_type ne 0, nwanted)
    nAllFiles = ulong64(total(parttotal[good]))  ; over all files

    ; init & alloc
    new_data = reform(new_data, rank, nData) ; consistent treatmnt of all ranks
    new_npart = ulong64(new_npart)

    if not keyword_set(all_data) then begin ; 1st file
        all_data_empty = 1
        all_data = make_array(rank, nAllFiles, type=type, val=0)
        all_npart = replicate(ulong64(0), 6) 
    end

    ; this is officially idx hell
    i_this = ulong64(0)                 ; i -> data

    j_this = ulong64(all_npart[0])      ; j -> all_data
    j_last = ulong64(total(all_npart))

    for type=0,5 do begin

        if want_type[type] ne 0 and new_npart[type] ne 0 then begin

            ; move old stuff out of the way
            if not keyword_set(all_data_empty) and nwanted gt 1 then begin 
                src_beg = ulong64(j_this)
                src_end = ulong64(j_last - 1)

                dst_beg = ulong64(j_this + new_npart[type])
                dst_end = ulong64(dst_beg + src_end - src_beg)
                
                for r=0, rank-1 do $
                    all_data[r, dst_beg:dst_end] = all_data[r, src_beg:src_end]
            end

            ; copy from read data
            dst_beg = j_this
            dst_end = dst_beg + new_npart[type] - 1

            src_beg = i_this
            src_end = src_beg + new_npart[type] - 1

            for r=0, rank-1 do $
                all_data[r, dst_beg:dst_end] = new_data[r, src_beg:src_end]

            ; update pointers
            all_npart[type] += ulong64(new_npart[type])
            j_this = ulong64(total(all_npart[0:type]))
            j_last = ulong64(total(all_npart))
        end
        
        ; update data pointer
        i_this += new_npart[type]

    end

    return
end

function read_gadget_snapshot, fname, label, header=header $
    , parttype=parttype, debug=debug
    
    close,/all

    if not keyword_set(fname) then begin
        print, "ReadGadgetSnap, fname, block, head=head, part=part"
        print, "    fname   : input file name"
        print, "    block   : block name, 4 chars long"
        print, "    head    : header"
        print, "    part    : specify type"
        print, "    debug   : print debug info / contents"
        return, -1
    end
    
    ; if you wanna make this a bitmask, go ahead
    want_type = [0,0,0,0,0,0]
    if n_elements(parttype) eq 0 then want_type[*] = 1 $
                                 else want_type[parttype] = 1

    ; make sure label is 4 chars long
    label = string((byte(label+"   "))[0:3]) 

    nfiles = get_nfiles(fname)

    if nfiles le 0 then begin
        print, 'Error: File "'+fname+'" not found !'
        return, -1
    end

    npart_read = make_array(6, /ulong, val=0)

    for file = 0,nFiles-1 do begin
        
        fin = fname
        if nFiles gt 1 then $
            fin += '.'+strn(file) 

        openr, fp, fin, /f77, /get_lun, $
            swap_endian=test_endianess(fin, debug=debug)
        
        header = read_gadget_header(fp)
        if label eq "HEAD" then $
            return, header

        parttotal = ulong64(header.parttotal)
        if nFiles eq 1 then $   ; read only one of many
            parttotal = ulong64(header.npart)

        blockinfo = get_block_info(label, header.npart, debug=debug)

        want_type *= blockinfo.ptype
        if total(want_type) eq 0 then begin ; check consistency
            print, 'Your types do not carry block '+label
            stop
        end

        if keyword_set(debug) then begin
            print, '----------------------------------' 
            print, 'File :', fin
            print, "npart:", header.npart
            print, "nAll :", header.parttotal
            print, "want :", want_type
        end

        while not eof(fp) do begin
            this_label = get_next_block_label(fp, blocksize=blocksize)

            if keyword_set(debug) then $
                print, this_label+" <-> "+label+', '+strn(blocksize)
            
            if this_label eq label then begin
                data = read_block(fp, blockinfo) 
                break    
            end 
            
            skip_block, fp, blocksize

        end

        free_lun, fp, /force
     
        ; catch special cases
        if label eq "MASS" and total(header.massarr) ne 0 then $
            data = generate_mass_from_header(header)

        if not keyword_set(data) then begin
            print, "Block "+label+' not found !'
            return, 0
        end

        ; rearrange particles
        sort_particles, data, header.npart, all_data=all_data, $
            all_npart=npart_read, want_type=want_type, parttotal=parttotal

        data = 0
    end
    
    heap_gc

    return, reform(all_data)    ; remove extra dim in rank=1 case
end



