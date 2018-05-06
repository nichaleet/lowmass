pro make_target_list
  spec = mrdfits('/scr2/nichal/workspace2/sps_fit/data/lowmass_lowsn_all/sps_fit.fits.gz',1)
  nobj = n_elements(spec)
  if file_test('/scr2/nichal/workspace2/lowmass/ana/lowmass_lowsnall_nichaleet_match.fits') eq 0 then begin
      cat = mrdfits('/scr2/nichal/workspace2/lowmass/casjobs/lowmass_lowsn_nichaleet.fit',1)
      newcat = replicate(cat(0),nobj)
      for i=0,nobj-1 do begin
         loc = where(cat.plate eq spec(i).plate and cat.fiberid eq spec(i).fiber and cat.mjd eq spec(i).mjd,nloc)
         if nloc ne 1 then stop
         newcat(i) = cat(loc)
      endfor
      mwrfits,newcat,'/scr2/nichal/workspace2/lowmass/ana/lowmass_lowsnall_nichaleet_match.fits',/create,/silent
   endif
   cat = mrdfits('/scr2/nichal/workspace2/lowmass/ana/lowmass_lowsnall_nichaleet_match.fits',1)
   restore,'/scr2/nichal/workspace2/lowmass/ana/lowmass_lowsnall_Ha_ew.sav'
   good = where(ha_ew_arr gt 0,cgood)
   plot,cat(good).ra/15.,cat(good).dec,psym=1
;lowleft 8 0 0 +0 0 0
;lowright 15 0 0 +0 0 0
;upleft 8 0 0 +60 0 0
;upright 15 0 0 +60 0 0
;most 13 0 0 +30 0 0
   good = where(ha_ew_arr gt 0 and cat.ra gt 90 and cat.ra lt 300 and cat.snmedian lt 15.,cgood)
   cat = cat(good)

   cat2 = mrdfits('/scr2/nichal/workspace2/lowmass/ana/lowmass_nichaleet_match.fits',1)
   cat2ra = mrdfits('/scr2/nichal/workspace2/lowmass/casjobs/lowmass_radec_nichaleet.fit',1)
   match = intarr(n_elements(cat2))
   for i=0,n_elements(Cat2)-1 do match(i) = where(cat2ra.plate eq cat2(i).plate and $
        cat2ra.fiberid eq cat2(i).fiberid and cat2ra.mjd eq cat2(i).mjd)
   cat2ra = cat2ra(match)
   if total(cat2ra.fiberid-cat2.fiberid) ne 0 then stop,'ohoh'
   restore,'/scr2/nichal/workspace2/lowmass/ana/lowmass_Ha_ew.sav'
   good = where(ha_ew_arr gt 0 and cat2ra.ra gt 90 and cat2ra.ra lt 300,cgood)
   cat2 = cat2(good)
   cat2ra = cat2ra(good)

   plot,cat(good).ra/15.,cat(good).dec,psym=1



   plot,cat.ra/15.,cat.dec,psym=1
   ;select 10 for each 0.1 logmass
   for i=0,11 do begin
      goodmass = where(cat.logmass gt 8+0.1*i and cat.logmass le 8.1+0.1*i,cgoodmass)
      if cgoodmass gt 10 then goodmass=goodmass[0:9]
      if i eq 0 then cattab = cat(goodmass) else cattab=[cattab,cat(goodmass)]
   endfor

 ;  for i=12,16 do begin
 ;     goodmass = where(cat2.logmass gt 8+0.1*i and cat2.logmass le 8.1+0.1*i,cgoodmass)
 ;     if cgoodmass gt 10 then goodmass=goodmass[0:9]
 ;     if i eq 12 then loc = goodmass else loc = [loc,goodmass]
 ;  endfor

   make_catalog_proposal,[cattab.ra,cat2ra(loc).ra],[cattab.dec,cat2ra(loc).dec],[cattab.logmass,cat2(loc).logmass],[cattab.PETROMAG_G,cat2(loc).petromag_g],[cattab.snmedian,cat2(loc).snmedian]

   cat = cattab
   set_plot,'ps'
   !p.font = 0
   psname = 'plots/radec.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plot,cat.ra/15.,cat.dec,xtitle='RA',ytitle='DEC',psym=cgsymcat(46)
   device,/close

   print, 'total ',nobj, 'with ',cgood, 'Ha'
   psname = 'plots/target_magnitude_mass_lowsn.eps'
   device, filename = psname,xsize = 12,ysize = 15, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   !p.multi=[0,1,2]
   vsym,4,/fill,rot=45
   plot,cat.petromag_g,cat.logmass,psym=8,xtitle='g mag',ytitle='log mass',yrange=[8.0,9.],/nodata,xrange=[15,19]
   cgerrplot,cat.petromag_g,cat.logmass_min,cat.logmass_max
   oplot,cat.petromag_g,cat.logmass,psym=8
   plothist,cat.petromag_g,bin=0.5,xrange=[15,19]
   device,/close

   psname = 'plots/target_mass_sn_lowsn.eps'
   !p.multi = [0,1,1]
   device, filename = psname,xsize = 15,ysize = 10, $
      xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   plot,cat.logmass,cat.snmedian,psym=8,xtitle='logmass',ytitle='sn'
   device,/close

   stop

end
