# Snapshot reading

# define the I/O blocks for the gadget format 2 file format
# and all known blocks with their data metrics

type IOBlock
	label::String		# 4 Chars exactly !
	name::String		#
	rank::Integer		# vector or scalar ?
	pType::UInt8 		# bitmask !
	elType::DataType	#
end

const IOBlocks = [	
				# All particles
				IOBlock("POS ", "Positions       ", 3, 0x7, Float32),
				IOBlock("VEL ", "Velocities      ", 3, 0x7, Float32),
			 	IOBlock("ID  ", "Identifiers     ", 1, 0x7, UInt32 ), 
				IOBlock("MASS", "Masses          ", 1, 0x7, Float32),
	 	
				# SPH only 
				IOBlock("RHO ", "Density         ", 1, 0x1, Float32),
				IOBlock("HSML", "Smoothing Length", 1, 0x1, Float32),
				IOBlock("U   ", "Internal Energy ", 1, 0x1, Float32), 
				IOBlock("BFLD", "Magnetic Field  ", 3, 0x1, Float32), 
				IOBlock("MACH", "Mach Number     ", 1, 0x1, Float32), 
				IOBlock("SHSP", "Shock Speed     ", 1, 0x1, Float32), 

				# Add above
				IOBlock("LAST", "Unknown Block   ", 1, 0x1, Float32) 
			   	]

type Header 
	npart::Array{UInt32}
	masses::Array{Float64}
	time::Float64
	redshift::Float64
	flag_sfr::Int32
	flag_feedback::Int32
	npartTotal::Array{UInt32}
	flag_cooling::Int32
	nFiles::UInt32
	boxsize::Float64
	omega0::Float64
	omegaL::Float64
	hbpar::Float64
	flag_comoving::Int32
	flag_double::Int32
	flag_periodic::Int32
	la::Array{UInt8} # == 256 Bytes

	function Header()

		npart = zeros(UInt32, 6)
		masses = zeros(Float64, 6)
		fill256 = zeros(UInt8, 84)

		return new(npart, masses, 0,0,0,0, npart, 0,0,0,0,0,0,0,0,0, fill256)
	end
end

function ReadSnap(fname::AbstractString, label::String; pType=0x7, debug=false)
	
	label = (label*"    ")[1:4] # 4 bytes

	nFiles = GetNFiles(fname)

	@assert(nFiles > 0, "\n\n    File not found $fname, $(fname).0\n")
	@assert(nFiles==1, "\n\n    Multiple file reading not implemented yet !\n")

	fd = open(fname, "r")
	
	head = ReadHeader(fd)

	if label == "HEAD"
		return head
	end

	blocksize = FindBlock(fd, label; debug=debug) # seek fd to block

	if blocksize == 0 && label == "MASS"
	
		data = MakeMassesFromHeader(head)

	else
		@assert(blocksize > 0, 
		  		"\n\n    Block <$label> not found in file '$fname' \n")

		data = ReadBlock(fd, label, blocksize; debug=debug)

		data = constrain!(data, head.npart, pType; debug=debug)

	end

	return data
end

function constrain!(data, npart, pType=0x07; debug=false)
	
	if pType == 0x07
		return data
	end

	cum = cumsum(npart)
  
	last = cum[pType+1]
	first = (1.+[0;cum])[pType+1]

	if debug == true
		println("Constraining block to $first : $last ")
	end

	if (size(data))[1] == 1
		return data[first:last]
	else
		return data[:, first:last]
	end

end

function ReadHead(fname::AbstractString; debug=false)

	return ReadSnap(fname, "HEAD"; debug=debug)
end


function ReadHeader(fd)

	seekstart(fd)

	label, nBytes = ReadFormat2Head(fd)

	@assert(label=="HEAD", "File does not begin with HEAD but $label")

	f77Record = read(fd, UInt32)

	@assert(f77Record == 256, "Error in F77 Record: $f77Record")
	
	head = Header()

	head.npart = read(fd, UInt32, 6)
	head.masses = read(fd, Float64, 6)
	head.time = read(fd, Float64)
	head.redshift = read(fd, Float64)
	head.flag_sfr = read(fd, Int32)
	head.flag_feedback = read(fd, Int32)
	head.npartTotal = read(fd, UInt32, 6)
	head.flag_cooling = read(fd, Int32)
	head.nFiles = read(fd, UInt32)
	head.boxsize = read(fd, Float64)
	head.omega0 = read(fd, Float64)
	head.omegaL = read(fd, Float64)
	head.hbpar = read(fd, Float64)
	head.flag_comoving = read(fd, Int32)
	head.flag_double = read(fd, Int32)
	head.flag_periodic = read(fd, Int32)
	head.la = read(fd, UInt8, 84)
	
	f77Record = read(fd, UInt32)

	@assert(f77Record == 256, "Error in F77 Record")

	return head
end

function FindBlock(fd, label::String; debug=false)

	seekstart(fd)

	blocksize = UInt64(0)

	while !eof(fd)

		blocklabel, skipsize = ReadFormat2Head(fd)

		blocksize = skipsize - 8 # account for F77 record

		if debug == true
			println("$label <=> $blocklabel, $(blocksize)")
		end

		if blocklabel == label

			break
		end

		skip(fd, skipsize)
	end
	
	return blocksize
end

function ReadBlock(fd::IO, label::String, blocksize;debug=debug)
	
	f77Record = read(fd, UInt32)

	@assert(f77Record == blocksize, 
		 "Fmt2 head size $blocksize does not match F77 record $f77Record")

	block = FindIOBlockFromLabel(label, debug=debug)

	npart = UInt64(blocksize / block.rank / sizeof(block.elType))

	if rank == 1
		data = zeros(block.elType, npart)
	else 
		data = zeros(block.elType, block.rank, npart)
	end

	read!(fd, data)

	f77Record = read(fd, UInt32)
		
	return data
end

function FindIOBlockFromLabel(label::String; debug=false)
	
	block = IOBlocks[1]

	for block in IOBlocks

		if label == block.label
			break
		end

	end

	if (debug == true) || (block.label == "LAST")

		println("Assuming '$(block.name)':")
		println("   label - $(block.label)")  
	 	println("   rank  - $(block.rank)")
		println("   ptype - $(block.pType)")
		println("   dtype - $(block.elType) ")
	end


	return block
end


function ReadFormat2Head(fd::IO)

	f77Record = read(fd, UInt32)
	@assert(f77Record == 8, "Error in F77 Record" )
	
	label = String(read(fd, UInt8, 4))
	nBytes = read(fd, UInt32) # size of block include f77 record (2*4bytes)

	f77Record = read(fd, UInt32)
	
	return label, nBytes
end

function WriteFormat2Head(fd::IO, label::String, blocksize::UInt32)

	@assert(length(label) == 4, "Label not 4 bytes long: $label")

	skipsize = UInt32(blocksize+8) 	# add F77 header size (2*4bytes)

	eight = UInt32(8)

	write(fd, eight)
	write(fd, label) 				# 4*char = 32 bits = 4 bytes
	write(fd, skipsize)				# 1*UInt32 = 32 bits = 4 bytes
	write(fd, eight)

end

function GetNFiles(fname::String)

	if isfile(fname)

		return 1

	elseif isfile(fname*".0")

		i = 0 # start with 0

		while isfile(fname*"."*string(i))
			i += 1
		end

		return i
	end

	return 0
end

function MakeMassesFromHeader(head::Header)

	nPart = sum(head.npart)

	data = Array{Float32}(nPart)

	run = 1

	for i=1:6

		for j=1:head.npart[i]
			
			data[run] = head.masses[i] 
	
			run += 1

		end

	end

	return data
end

function WriteHead(fname::AbstractString, head::Header; replace=false, debug=false)

	if debug == true
		print("Writing HEAD to file '$fname' ")
	end

	if (replace==true) && isfile(fname)
		rm(fname)
	end
	
	@assert(!isfile(fname) || (replace==true), 
		 	"\n   File '$fname' exists already and replace==false !\n\n" ) 

	fd = open(fname, "a+")

	seekstart(fd)

	WriteFormat2Head(fd, "HEAD", UInt32(256))

	f77record = UInt32(256)

	write(fd, f77record)

	write(fd, head.npart)
	write(fd, head.masses)
	write(fd, head.time)
	write(fd, head.redshift)
	write(fd, head.flag_sfr)
	write(fd, head.flag_feedback)
	write(fd, head.npartTotal)
	write(fd, head.flag_cooling)
	write(fd, head.nFiles)
	write(fd, head.boxsize)
	write(fd, head.omega0)
	write(fd, head.omegaL)
	write(fd, head.hbpar)
	write(fd, head.flag_comoving)
	write(fd, head.flag_double)
	write(fd, head.flag_periodic)
	write(fd, head.la) # == 256 bytes

	write(fd, f77record)

	close(fd)
	
	if debug == true
		println("done")
	end
end

function AddBlock(fname::AbstractString, label::String, data; debug=false)

	dType = typeof(data)
	blocksize = UInt32(sizeof(data))
	rank = (size(data))[1]
	npart = (size(data))[2]

	label *= "    "
	label = label[1:4]

	if debug == true
		println("Adding block $label to file $fname")
		println("   rank  = $rank")
		println("   npart = $npart")
		println("   type  = $dType")
	end

	fd = open(fname, "a+")

	seekstart(fd)

	while !eof(fd) # skip manually to avoid double blocks

		fLabel, skipsize = ReadFormat2Head(fd)

		@assert(fLabel != label, "Label $label already exists in file $fname")

		if debug == true
			println("   Skipping fwd: block $fLabel <-> $skipsize")
		end

		skip(fd, skipsize)
	end

	# hopefully at the end now

	WriteFormat2Head(fd, label, blocksize)

	f77record = UInt32(blocksize)

	write(fd, f77record)

	write(fd, data)
	
	write(fd, f77record)

	close(fd)

end
