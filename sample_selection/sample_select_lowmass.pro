function match_spec_and_casjob,specdir,casjobfile
   casstr = mrdfits(casjobfile,1)
   specfiles = file_search(specdir+'spec-*.fits',count=nspec)

   outstr = create_struct(casstr[0],'specfile','fileloc','haew',0.,'haew_err',0.,$
                          'kcorrectmass',0.,'kcorrectmass_dered',0.,'itime',0.,'SN_itime',0.)
   outstrarr = replicate(outstr,nspec) 
   for i=0,nspec-1 do begin
      specinfo = mrdfits(specfiles[i],2,/silent)
      matchloc = where(casstr.plate eq specinfo.plate and casstr.fiberid eq specinfo.fiberid $
                       and casstr.mjd eq specinfo.mjd,cmatchloc)
      if cmatchloc eq 1 then begin
         casstrnow = casstr[matchloc]
         struct_assign,casstrnow,outstr
         outstrarr[i] = outstr
         outstrarr[i].specfile = specfiles[i]
      endif else stop,'no matched catalog found'          
   endfor   
   return, outstrarr
end

pro sample_select_lowmass,redomatch=redomatch,redoha=redoha,redokcorrect=redokcorrect,recalitime=recalitime
;this is the master run for sample selection
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;create a master files of all spectra from casjob with added tags of Ha, spectrum name
   matchedfits = 'lowmass_allsn_match_spec_casjob.fits'
   if file_test(matchedfits) eq 0 or keyword_set(redomatch) then begin 
      casjobdir = '/scr2/nichal/workspace_lowmass/casjobs/'
      matchedstr = match_spec_and_casjob(casjobdir+'spec_lowmass_allsn/',$
                   casjobdir+'lowmass_allsn_nichaleet.fit')
      mwrfits, matchedstr, matchedfits,/create,/silent
   endif
;calculate Ha EW
   if keyword_set(redoha) then begin
      mainstr = mrdfits(matchedfits,1)
      measurehaew,mainstr,save_every=10,savename=matchedfits
      mwrfits,mainstr, matchedfits,/create,/silent
   endif
;calculate kcorrect mass
   if keyword_set(redokcorrect) then begin
      mainstr = mrdfits(matchedfits,1)
      get_kcorrectmass,mainstr
      mwrfits,mainstr,matchedfits,/create,/silent
   endif
;calculate integration time
   if keyword_set(recalitime) then begin
      mainstr = mrdfits(matchedfits,1)
      sn = 30.
      mainstr.itime = calitime_dbsp(mainstr.petromag_g,sn,4700.,600.) 
      mainstr.sn_itime = sn
      mwrfits,mainstr,matchedfits,/create,/silent
   endif
   mainstr = mrdfits(matchedfits,1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   
;select only spectrum with Ha EW less than 0 
   passive = where(mainstr.haew lt 0. and finite(mainstr.haew_err) and mainstr.haew gt -20.,cpassive)
   strnow = mainstr(passive)
   ;plothist,strnow.logmass

;high priority select all with log(mass) < 8.7 with RA good for both april and may obs
   lowmass_group = where(strnow.logmass lt 8.7 and strnow.ra gt 130.,clowmass_group)
   strlow = strnow(lowmass_group)
   strlow = strlow(sort(strlow.ra))
   ;sort by ra then writeout an output
   index = indgen(clowmass_group)+1
   writecol,'lowmass_targets.txt',index,strlow.ra,strlow.dec,strlow.z,strlow.logmass,$
             strlow.petromag_g,strlow.petromag_r,strlow.snmedian,strlow.haew,strlow.itime/60.,$
             fmt='(I03,x,f11.7,x,f11.7,x,f6.4,x,f5.3,x,f5.2,x,f5.2,x,f5.2,x,f6.2,x,f5.2)'
   writecol,'lowmass_objvisibility.txt',index,strlow.ra,strlow.dec,fmt='(I03,x,f11.7,x,f11.7)'
   mwrfits,strlow,'lowmass_targets.fits',/create,/silent
stop
end
