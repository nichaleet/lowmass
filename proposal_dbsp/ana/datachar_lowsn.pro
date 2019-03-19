pro datachar_lowsn
  spec = mrdfits('/scr2/nichal/workspace2/sps_fit/data/lowmass_lowsn/sps_fit.fits.gz',1)
  nobj = n_elements(spec)
  if file_test('/scr2/nichal/workspace2/lowmass/ana/lowmass_lowsn_nichaleet_match.fits') eq 0 then begin
      cat = mrdfits('/scr2/nichal/workspace2/lowmass/casjobs/lowmass_lowsn_nichaleet.fit',1)
      newcat = replicate(cat(0),nobj)
      for i=0,nobj-1 do begin
         loc = where(cat.plate eq spec(i).plate and cat.fiberid eq spec(i).fiber and cat.mjd eq spec(i).mjd,nloc)
         if nloc ne 1 then stop
         newcat(i) = cat(loc)
      endfor
      mwrfits,newcat,'/scr2/nichal/workspace2/lowmass/ana/lowmass_lowsn_nichaleet_match.fits',/create,/silent
   endif
   cat = mrdfits('/scr2/nichal/workspace2/lowmass/ana/lowmass_lowsn_nichaleet_match.fits',1)
   restore,'lowmass_lowsn_Ha_ew.sav'

   good = where(ha_ew_arr gt 0,cgood)
   spec.good = 0
   spec(good).good = 1
   spec.logmstar = cat.logmass
   mwrfits,spec,'/scr2/nichal/workspace2/sps_fit/data/lowmass_lowsn/sps_fit.fits.gz',/create,/silent

   print, 'total ',nobj, 'with ',cgood, 'mass'

   set_plot,'ps'
   !p.font = 0
   psname = 'plots/mass_hist_lowsn.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plothist,cat.logmass,bin=0.1
   plothist,cat(good).logmass,color=fsc_Color('red'),/overplot,bin=0.1
   device,/close

   psname = 'plots/magnitude_mass_lowsn.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   vsym,4,/fill,rot=45
   plot,cat.petromag_g,cat.logmass,psym=8,xtitle='g mag',ytitle='log mass',yrange=[8.0,9.3]
   cgerrplot,cat.petromag_g,cat.logmass_min,cat.logmass_max
   oplot,cat(good).petromag_g,cat(good).logmass,psym=8,color=fsc_color('maroon')
   device,/close

   psname = 'plots/mass_sn_lowsn.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plot,cat.logmass,cat.snmedian,psym=8,xtitle='logmass',ytitle='sn'
   oplot,cat(good).logmass,cat.snmedian,psym=8,color=fsc_color('maroon')
   device,/close

   stop

end
