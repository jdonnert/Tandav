pro generate_fileio_data
	
	gadget = obj_new('GADGETCODEOBJECT')

	outfile = "./data/IC_fileio"

	nFiles = 12

	npart = 100000

	for file=0, nFiles-1 do begin
		
		head = gadget.makeHead()

		head.npart[1] = npart
		head.parttotal = head.npart * nFiles

		head.num_files = nFiles
		head.redshift = 1
		head.time = 1 / (1+head.redshift)

		fout = outfile+"."+strn(file, len=3, padc='0')

		gadget.writeHead, fout, head
		
		data = make_array(3, npart)
		data[0,*] = findgen(npart)
		data[1,*] = findgen(npart)
		data[2,*] = findgen(npart)
		gadget.addBlock, fout, float(data), 'POS'

		data = make_array(3, npart)
		data[0,*] = findgen(npart)
		data[1,*] = findgen(npart)
		data[2,*] = findgen(npart)
		gadget.addBlock, fout, float(data), 'VEL'

		data = 1+lindgen(npart) + file*npart
		gadget.addBlock, fout, ulong(data), 'ID'

		data = make_array( npart)
		gadget.addBlock, fout, float(data), 'MASS'
	end

	return
end
