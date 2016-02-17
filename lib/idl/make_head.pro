; make a tandav header
function make_head, npart=npart
	head = {npart		    :	lonarr(6),	$
			massarr			:	dblarr(6),	$
			time			: 	double(1),	$
			redshift		:	double(1),	$
			flag_sfr		:	long(0),	$
			flag_feedback	:	long(0), 	$
			parttotal		: 	lonarr(6),	$
			flag_cooling	:	long(0),	$
			num_files		:	long(1),	$
			boxsize			:	double(1),	$
			omega0			:	double(0.3),$
			omegalambda		:	double(0.7),$
			hubbleparam		:	double(0.7),$
			la				:	make_array(96,/byte) }

	return, head
end

