pro sample_select_highermass,redoaddsdsshaew=redoaddsdsshaew
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;use the matchedfits from sample_select_lowmass
   matchedfits = 'lowmass_allsn_match_spec_casjob.fits'
   if keyword_set(redoaddsdsshaew) then begin
      mainstr = mrdfits(matchedfits,1)
      nobj = n_elements(mainstr)
      newstr = []
      for i=0,nobj-1 do begin
          strnow = mainstr[i]
          struct_add_field,strnow,'haew_sdss',0.
          struct_add_field,strnow,'haewerr_sdss',0.
          struct_add_field,strnow,'sdss_exptime',0.
          sdsslines  = mrdfits(strnow.specfile,3,/silent)
          loc = where(strmatch(sdsslines.linename,'*H_alpha*'),cloc)
          if cloc eq 1 then begin
             strnow.haew_sdss = sdsslines[loc].lineew
             strnow.haewerr_sdss = sdsslines[loc].lineew_err
          endif
          fits_info,strnow.specfile,/silent,n_Ext=n_ext
          totexptime=0.           
          for j=4,n_ext do begin
             aa=mrdfits(strnow.specfile,j,hdr,/silent)
             loc = where(strmatch(hdr,'*EXPTIME*'),cloc)
             if cloc eq 1 then begin
               aa=strsplit(hdr(loc),/extract)
               totexptime += float(aa[2])
             endif
          endfor
          strnow.sdss_Exptime = totexptime
          newstr = [newstr,strnow]
      endfor
      mwrfits,newstr, matchedfits,/create,/silent
   endif
   mainstr = mrdfits(matchedfits,1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   
;select only spectrum with Ha EW less than 0 
   passive = where(mainstr.haew lt 0. and finite(mainstr.haew_err) and mainstr.haew gt -20. and mainstr.haew_sdss lt 0.,cpassive)
   strnow = mainstr(passive)
;   p = plot(strnow.logmass,strnow.snmedian,xtitle='mass',ytitle='SDSS SN',symbol='star',linestyle=6)
;   p=plot(!x.crange,[25,25],/overplot,color='red')

;   good = where(strnow.snmedian gt 25.,cgood)
;   hist = histogram(strnow(good).logmass,binsize=0.1,locations=xbin)
;   phist = plot(xbin,hist,/histogram,xtitle='mass')

;lower priority select all with 8.7 <= log(mass) < 8.9 with RA good for both april and may obs
;need 2 of log(mass) between 8.8 and 8.9 to make sample of 5 with SN>25 in this mass bin
;need 3 of log(mass) between 8.7 and 8.8

   highermass_group = where(strnow.logmass lt 8.9 and strnow.logmass ge 8.7 and strnow.ra gt 130. and strnow.snmedian gt 10 and strnow.snmedian lt 25.,chighermass_group)
   strhigh = strnow(highermass_group)
   strhigh = strhigh(sort(strhigh.ra))
   ;sort by ra then writeout an output
   index = indgen(chighermass_group)+74
   writecol,'highermass_targets.txt',index,strhigh.ra,strhigh.dec,strhigh.z,strhigh.logmass,$
             strhigh.petromag_g,strhigh.petromag_r,strhigh.snmedian,strhigh.haew,strhigh.itime/60.,$
             fmt='(I03,x,f11.7,x,f11.7,x,f6.4,x,f5.3,x,f5.2,x,f5.2,x,f5.2,x,f6.2,x,f5.2)'
   writecol,'highermass_objvisibility.txt',index,strhigh.ra,strhigh.dec,fmt='(I03,x,f11.7,x,f11.7)'
   mwrfits,strhigh,'highermass_targets.fits',/create,/silent
stop
end
