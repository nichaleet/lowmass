pro test
   str = mrdfits('highermass_targets.fits',1) 
   file = str[0].specfile
   fits_info,file,/silent,n_Ext=n_ext 
   for i=4,n_ext do begin
      aa=mrdfits(file,i,hdr,/silent)
      print, i
      print, hdr[[21,25,33]]
      stop
   endfor
stop
end
