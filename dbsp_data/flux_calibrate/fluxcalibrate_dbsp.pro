pro fluxcalibrate_dbsp,stdfile,obsfile
   ;output is calibration curve which is real_flambda = curve*obs_flambda
   ;file out is stdfile+.calib.fits

   stddir = '/scr2/nichal/workspace_lowmass/dbsp_data/flux_calibrate/standard_stars/'
   outfile = (strsplit(file_basename(obsfile),'.',/extract))[0]+'.calib.fits'
   ;read standard star values
   readcol,stddir+stdfile,stdwl,stdmag,stdsomething
   ;change abmag to flambda
   stdfnu = 10.^(-0.4*(stdmag+48.6))
   stdflambda = 3.e8/(stdwl*1.e-8)^2*stdfnu

   ;read observations
   flux = readfits(obsfile,hdr)
   cdelt1 = sxpar(hdr,'cdelt1')
   crval1 = sxpar(hdr,'crval1')
   lambda = findgen(n_elements(flux))*cdelt1+crval1

   ;interpolate standard to observations
   stdflambda = interpol(stdflambda,stdwl,lambda)
   
   curve = stdflambda/flux
   strout={standard_file:stddir+stdfile,fcalib:curve,lambda:lambda,hdr:hdr}
   mwrfits,strout,outfile,/create,/silent

end
   
