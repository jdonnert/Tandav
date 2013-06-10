pro write_head, fname, h, swap_endian=swap_endian
   bl = long(264)

   b4 = byte("HEAD    ")
   b4 = b4(0:3)

   openw, fd, fname, /f77, /get_lun, swap_endian=swap_endian

   writeu,fd, b4, bl
   
   writeu, fd, h.npart, h.massarr, h.time, h.redshift, h.flag_sfr, $
       h.flag_feedback, h.partTotal, h.flag_cooling, h.num_files, $
       h.BoxSize, h.Omega0, h.OmegaLambda, h.HubbleParam, h.la
   
   close, fd

   free_lun, fd
END   
