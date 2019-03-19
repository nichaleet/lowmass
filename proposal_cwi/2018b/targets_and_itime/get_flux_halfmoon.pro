pro get_flux_halfmoon
;take the sample from sample_selections found for dbsp
;calculate flux per arcsecond and area
   cat1lowmass = mrdfits('/scr2/nichal/workspace_lowmass/sample_selection/lowmass_targets.fits',1)
   cat1himass = mrdfits('/scr2/nichal/workspace_lowmass/sample_selection/highermass_targets.fits',1)
   cat1 = [cat1lowmass,cat1himass]

   redomatch = 0
   if redomatch then begin
      cat2 = mrdfits('petromag_info.fits',1);this was downloaded from casjob
      ngals = n_elements(cat1)
     
      loccat2 = lonarr(ngals)
      for i=0,ngals-1 do begin
         loc = where(cat1[i].plate eq cat2.plate and cat1[i].mjd eq cat2.mjd and $
                     cat1[i].fiberid eq cat2.fiberid,cloc)
         if cloc ne 1 then stop,'something wrong'
         loccat2[i] = loc
      endfor
      cat2 = cat2(loccat2)
      mwrfits,cat2,'petromag_info_matched_targets.fits',/create,/silent
   endif
   cat2 = mrdfits('petromag_info_matched_targets.fits',1)
   ngals = n_elements(cat2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;calculate flux and itime
redoitime = 1
if redoitime eq 1 then begin
   clight = 2.99792458d10 ;cgs
   lambda = 4770. ;Angstrom gband
   conv = 3.631d-6*1.d-23*clight*1.d8/lambda^2 ;change nanomaggie to erg/s/cm2/A
   ;nanomaggie; 1 nanomaggie= 3.631×10-6 Jy, 1 Jy = 1.d-23 cgs
   strtemp = {plate:0s,mjd:0L,fiberid:0s,petro50flux:0.,petro50fluxerr:0.,petro90flux:0.,$
          petro90fluxerr:0.,petro50itime:0.,petro90itime:0.,petro50itime_kcwi:0.,petro90itime_kcwi:0.,petrooutringitime:0.}
   strout = replicate(strtemp,ngals)
   strout.plate = cat2.plate
   strout.mjd = cat2.mjd
   strout.fiberid = cat2.fiberid
   for i=0,ngals-1 do begin
      petroflux = (cat2[i].petroflux_g)*conv ;erg/s/cm2/A
      petrofluxerr = 1.d/sqrt(cat2[i].petrofluxivar_g)*conv
      area50 = !const.pi*(cat2[i].PETROR50_G)^2
      area50err = !const.pi*2.*cat2[i].PETROR50_G*cat2[i].PETROR50ERR_G
      area90 = !const.pi*(cat2[i].PETROR90_G)^2
      area90err = !const.pi*2.*cat2[i].PETROR90_G*cat2[i].PETROR90ERR_G
      strout[i].petro50flux = petroflux*0.5/area50  ;erg/s/cm2/A/arcsec2
      strout[i].petro50fluxerr = strout[i].petro50flux*sqrt(petrofluxerr^2/petroflux^2+area50err^2/area50^2)
      strout[i].petro90flux = petroflux*0.9/area90  ;erg/s/cm2/A/arcsec2
      strout[i].petro90fluxerr = strout[i].petro90flux*sqrt(petrofluxerr^2/petroflux^2+area90err^2/area90^2)
      ;for 1.e-17 erg/s/cm2/A/arcsec^2 get s/n=5 per arcsec^2 in 1 hr
      strout[i].petro90itime = cwi_etc_nl_caltime(strout[i].petro90flux,25.,snr_spatial_bin=area90,/halfmoon)      
      strout[i].petro50itime = cwi_etc_nl_caltime(strout[i].petro50flux,25.,snr_spatial_bin=area50,/halfmoon)
      aveflux_outring = petroflux*0.4/(area90-area50)
      strout[i].petrooutringitime = cwi_etc_nl_caltime(aveflux_outring,15.,snr_spatial_bin=area90-area50,/halfmoon)
   endfor
   print, 'finish'
   mwrfits,strout,'itime_targets_halfmoon.fits',/create,/silent
endif 
end
