# Snapshot reading

# define the I/O blocks for the gadget format 2 file format
# and all known blocks with their data metrics

module TandavCodeSnapshotIO

export Header 
export ReadHead, ReadSnap

type IOBlock
	label::String		# 4 Chars exactly !
	name::String		#
	rank::Integer		# vector or scalar ?
	pType::UInt8 		# bitmask !
	elType::DataType	#
end

const IOBlocks = [	IOBlock("POS ", "Positions      ", 3, 0x7, Float32),
			 		IOBlock("VEL ", "Velocities     ", 3, 0x7, Float32),
				 	IOBlock("ID  ", "Identifiers    ", 1, 0x7, UInt32 ), 
				 	IOBlock("MASS", "Masses         ", 1, 0x7, Float32),
	 	
					IOBlock("RHO ", "Density        ", 1, 0x1, Float32),
				 	IOBlock("U   ", "Internal Energy", 1, 0x1, Float32), 
				 	IOBlock("BFLD", "Magnetic Field ", 3, 0x1, Float32), 
				 	IOBlock("MACH", "Mach Number    ", 1, 0x1, Float32), 
				 	IOBlock("LAST", "Unknown Block  ", 1, 0x1, Float32) 
			   	]

type Header 
	npart::Vector{UInt32}
	masses::Vector{Float64}
	time::Float64
	redshift::Float64
	flag_sfr::Int32
	flag_feedback::Int32
	npartTotal::Vector{UInt32}
	flag_cooling::Int32
	nFiles::UInt32
	boxsize::Float64
	omega0::Float64
	omegaL::Float64
	hbpar::Float64
	flag_comoving::Int32
	flag_double::Int32
	flag_periodic::Int32
	la::Vector{UInt8} # == 256 Bytes

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

	@assert(blocksize > 0, "\n\n    Block <$label> not found in file '$fname' \n")

	data = ReadBlock(fd, label, blocksize)

	return data
end

function ReadHead(fname::AbstractString; debug=false)

	return ReadSnap(fname, "HEAD"; debug=debug)
end


function ReadHeader(fd)

	seekstart(fd)

	label, nBytes = ReadFormat2Head(fd)

	@assert(label=="HEAD", "File does not begin with HEAD")

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

		if debug == true
			println("$label <=> $blocklabel, $skipsize")
		end

		if blocklabel == label

			blocksize = skipsize - 8 # account for F77 record

			break
		end

		skip(fd, skipsize)
	end
	
	return blocksize
end

function ReadBlock(fd::IO, label::String, blocksize)
	
	f77Record = read(fd, UInt32)

	@assert(f77Record == blocksize, "Fmt2 head size does not match F77 record")

	elType = Float32
	rank = 1

	for blk in IOBlocks

		if label == blk.label

			elType = blk.elType
			rank = blk.rank

			break
		end

		if blk.label == "LAST"
			println("Block $label not known. Assuming 1D SPH only.")
		end 
	end

	npart = UInt64(blocksize / rank / sizeof(elType))

	if rank == 1
		data = zeros(elType, npart)
	else 
		data = zeros(elType, rank, npart)
	end

	read!(fd, data)

	f77Record = read(fd, UInt32)

	return data
end


function ReadFormat2Head(fd::IO)

	f77Record = read(fd, UInt32)
	@assert(f77Record == 8, "Error in F77 Record" )
	
	label = String(read(fd, UInt8, 4))
	nBytes = read(fd, UInt32) # size of block include f77 record (2*4bytes)

	f77Record = read(fd, UInt32)
	
	return label, nBytes
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

function WriteHead(fname::AbstractString, head::Header)

end

function AddBlock(fname::AbstractString, label::String)


end
end # module
