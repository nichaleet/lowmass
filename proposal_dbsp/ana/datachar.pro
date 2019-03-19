pro datachar
   spec = mrdfits('/scr2/nichal/workspace2/sps_fit/data/lowmass/sps_fit.fits.gz',1)
   nobj = n_elements(spec)
   if file_test('/scr2/nichal/workspace2/lowmass/ana/lowmass_nichaleet_match.fits') eq 0 then begin
      cat = mrdfits('/scr2/nichal/workspace2/lowmass/casjobs/lowmass_nichaleet.fit',1)
      newcat = cat
      for i=0,nobj-1 do begin
         loc = where(cat.plate eq spec(i).plate and cat.fiberid eq spec(i).fiber and cat.mjd eq spec(i).mjd,nloc)
         if nloc ne 1 then stop
         newcat(i) = cat(loc)
      endfor
      mwrfits,newcat,'/scr2/nichal/workspace2/lowmass/ana/lowmass_nichaleet_match.fits',/create,/silent
   endif
   cat = mrdfits('/scr2/nichal/workspace2/lowmass/ana/lowmass_nichaleet_match.fits',1)
   restore,'lowmass_Ha_ew.sav'

   good = where(ha_ew_arr gt 0,cgood)
   spec.good = 0
   spec(good).good = 1
   spec.logmstar = cat.logmass
   mwrfits,spec,'/scr2/nichal/workspace2/sps_fit/data/lowmass/sps_fit.fits.gz',/create,/silent

   set_plot,'ps'
   !p.font = 0
   psname = 'plots/mass_hist.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plothist,cat.logmass,yrange=[0,60]
   plothist,cat(good).logmass,color=fsc_Color('red'),/overplot
   device,/close

   psname = 'plots/sn_sn.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plot,cat.snmedian,spec.sn,xtitle='sdss nominal sn',ytitle='measured sn',psym=1
   oplot,[0,50],[0,50],color=fsc_color('red')
   device,/close

   psname = 'plots/magnitude_mass.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   vsym,4,/fill,rot=45
   plot,cat.petromag_g,cat.logmass,psym=8,xtitle='g mag',ytitle='log mass',yrange=[8.2,9.5]
   cgerrplot,cat.petromag_g,cat.logmass_min,cat.logmass_max
   oplot,cat(good).petromag_g,cat(good).logmass,psym=8,color=fsc_color('red')
   device,/close

   psname = 'plots/mass_sn.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plot,cat.logmass,cat.snmedian,psym=4,xtitle='logmass',ytitle='sn'
   device,/close

   stop



end
