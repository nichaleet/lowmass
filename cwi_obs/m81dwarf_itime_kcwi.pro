pro m81dwarf_itime_kcwi

   clight = 2.99792458d10 ;cgs

   ;d99744
   sfb = 23.9;surface brightness within effective radius in mag per arcsec2 ;in r band
   lambda = 6400. ;angstrom
   snr_arr = [5,10,15,20,25,30]
   re = 18.3
   spatialbin = !pi*re^2
   half_spatialbin = !pi*(0.5*re)^2
   ;calculating flux
   sff  = 1.d8*clight*10.^((sfb+48.6)/(-2.5))/lambda^2 ; flambda in erg per (s cm2 angstrom arcsec2)
   print,'d09744 integration time'
   print,'surface brightness is ',sff, 'erg per (s cm2 angstrom arcsec2)'
   print, 'SNR in Re         Itime(Re)       Itime(half Re)'
   for i=0,n_elements(snr_arr)-1 do begin
       snr = snr_arr[i] 
       itime = kcwi_etc_nl_caltime(sff,snr,snr_spatial_bin=spatialbin)    
       itime2 = kcwi_etc_nl_caltime(sff,snr,snr_spatial_bin=half_spatialbin)    
       print, snr, itime,itime2
   endfor

   ;K61
   sfb = 24.7 ;within Re
   lambda = 6400. ;
   sfb_central = 23.76 ;central SB in V
   lambda_central = 5500. 
   re = 38.4
   spatialbin = !pi*re^2
   half_spatialbin = !pi*(0.5*re)^2
   sff  = 1.d8*clight*10.^((sfb+48.6)/(-2.5))/lambda^2 ; flambda in erg per (s cm2 angstrom arcsec2)
   sff_central  = 1.d8*clight*10.^((sfb_central+48.6)/(-2.5))/lambda_central^2 
   print,'KDG61 integration time'
   print,'surface brightness is ',sff, 'erg per (s cm2 angstrom arcsec2)'
   print,'central surface brightness is ',sff_central, 'erg per (s cm2 angstrom arcsec2)'
   print, 'SNR in Re         Itime(Re)       Itime(half Re)    Itime(central)'
   for i=0,n_elements(snr_arr)-1 do begin
       snr = snr_arr[i]
       itime = kcwi_etc_nl_caltime(sff,snr,snr_spatial_bin=spatialbin)
       itime2 = kcwi_etc_nl_caltime(sff,snr,snr_spatial_bin=half_spatialbin)
       itime3 = kcwi_etc_nl_caltime(sff_central,snr,snr_spatial_bin=half_spatialbin)
       print, snr, itime,itime2,itime3
   endfor

   ;K64
   sfb = 23.8 ;within Re
   lambda = 6400. ;
   sfb_central = 22.89 ;central SB in V
   lambda_central = 5500.
   re = 22.1
   spatialbin = !pi*re^2
   half_spatialbin = !pi*(0.5*re)^2
   sff  = 1.d8*clight*10.^((sfb+48.6)/(-2.5))/lambda^2 ; flambda in erg per (s cm2 angstrom arcsec2)
   sff_central  = 1.d8*clight*10.^((sfb_central+48.6)/(-2.5))/lambda_central^2
   print,'KDG64 integration time'
   print,'surface brightness is ',sff, 'erg per (s cm2 angstrom arcsec2)'
   print,'central surface brightness is ',sff_central, 'erg per (s cm2 angstrom arcsec2)'
   print, 'SNR in Re         Itime(Re)       Itime(half Re)    Itime(central)'
   for i=0,n_elements(snr_arr)-1 do begin
       snr = snr_arr[i]
       itime = kcwi_etc_nl_caltime(sff,snr,snr_spatial_bin=spatialbin)
       itime2 = kcwi_etc_nl_caltime(sff,snr,snr_spatial_bin=half_spatialbin)
       itime3 = kcwi_etc_nl_caltime(sff_central,snr,snr_spatial_bin=half_spatialbin)
       print, snr, itime,itime2,itime3
   endfor

  ;K63
   sfb = 24.1;surface brightness within effective radius in mag per arcsec2 ;in r band
   lambda = 6400. ;angstrom
   snr_arr = [5,10,15,20,25,30]
   re = 24.9
   spatialbin = !pi*re^2
   half_spatialbin = !pi*(0.5*re)^2
   ;calculating flux
   sff  = 1.d8*clight*10.^((sfb+48.6)/(-2.5))/lambda^2 ; flambda in erg per (s cm2 angstrom arcsec2)
   print,'KDG63 integration time'
   print,'surface brightness is ',sff, 'erg per (s cm2 angstrom arcsec2)'
   print, 'SNR in Re         Itime(Re)       Itime(half Re)'
   for i=0,n_elements(snr_arr)-1 do begin
       snr = snr_arr[i]
       itime = kcwi_etc_nl_caltime(sff,snr,snr_spatial_bin=spatialbin)
       itime2 = kcwi_etc_nl_caltime(sff,snr,snr_spatial_bin=half_spatialbin)
       print, snr, itime,itime2
   endfor

end
