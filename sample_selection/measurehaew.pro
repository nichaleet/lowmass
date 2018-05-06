pro measurehaew,casstr,save_every=save_every,savename=savename
   nspec = n_Elements(casstr)   
   set_plot,'x'
   for i=0,nspec-1 do begin
      print, 'doing '+sigfig(i,3)+'/'+sigfig(nspec,3)+casstr[i].specfile
      spec = mrdfits(casstr[i].specfile,1,/silent)
      zin = casstr[i].z
      flux = spec.flux
      lambda = 10.^spec.loglam
      ivar = spec.ivar
      casstr[i].haew = Ha_ew(flux,lambda,ivar,zin,widtherr,/nostop)
      casstr[i].haew_err = widtherr
      if keyword_set(save_every) then begin
         if i mod save_every eq 0 then begin
            mwrfits,casstr, savename,/create,/silent
            print, 'save ',i
         endif
      endif
   endfor
   stop
end
