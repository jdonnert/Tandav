pro generate_fileio_data
	
	gadget = obj_new('GADGETCODEOBJECT')

	outfile = "./data/IC_fileio"

	nFiles = 10

	npart = 4000000

	for file=0, nFiles-1 do begin
		
		head = gadget.makeHead()

		head.npart[0] = npart / 4.0
		head.npart[1] = npart / 4.0
		head.npart[2] = npart / 4.0
		head.npart[3] = npart / 4.0
		head.parttotal = head.npart * nFiles

		head.massarr = [0.23003, 0.23421, 0.423555, 0.44442, 0, 0]
		head.num_files = nFiles
		head.redshift = 1
		head.time = 1 / (1+head.redshift)

		fout = outfile+"."+strn(file)

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

		data = 1+lindgen(npart) + file*total(head.npart)
		gadget.addBlock, fout, ulong(data), 'ID'

		data = make_array( npart)
		gadget.addBlock, fout, float(data), 'MASS'
	end

	return
end
