pro quickspec,x,a,f,pder
common spscommon, sps, spsage
common dispersion, cdelt
   clight = 299792.458
   age = a[0]
   veldisp = a[1]
   redshift = a[2]
   lambdain = x/(redshift+1.)   
   nspsage = n_elements(spsage)
   ;interpolate for age
   wa = value_locate(spsage, age)
   if spsage[wa] eq age then begin
       ia = wa
       da = 1d
       na = 1
   endif else begin
       if wa eq nspsage-1 then wa -= 1
       ia = [wa, wa+1]
       da = abs(reverse(spsage[ia])-age)
       na = 2
   endelse
   spsspec = sps[0].spec*0.
   for j=0,na-1 do begin
       w = where(sps.agegyr eq spsage[ia[j]], c)
       if c eq 0 then message, 'Spectrum not found.'
       if j eq 0 then begin
          lambda = sps[w].lambda
          spsspeci = sps[w].spec
       endif else spsspeci = interpol(sps[w].spec, sps[w].lambda, lambda)
       spsspec += da[j]*spsspeci
   endfor
   ;interpolate to lambdain
   spsspec = interpol(spsspec,lambda,lambdain)

   ;smooth to velocity dispersion
   sigma = veldisp/clight*median(lambdain)/cdelt ;pixels
   spsspec = gauss_smooth(spsspec,sigma,/edge_truncate)
  
   ;continuum normalized
   contpara = poly_fit(x,spsspec,6,yfit=cont)
   spsspec = spsspec/cont
   f = spsspec
end

pro quicksn_dbsp,file,redshift
common spscommon, sps, spsage
common dispersion, cdelt
;quickly fit templates to reduced data and find signal to noise
;metallicity is fixed at -0.4 dex, only fit for age,redshift, and velocity dispersion
;used to calculate S/N while observing
   spec = readfits(file,hdr)
   cdelt = sxpar(hdr,'cdelt1')
   lambda = (findgen(n_elements(spec))-sxpar(hdr,'crpix1')+1.)*cdelt+sxpar(hdr,'crval1')
   ;continuum normalize
   ;polynomial degree 6 (same as choi14)
   contpara = poly_fit(lambda,spec,6,yfit=cont)
   plot,lambda,spec,xrange=minmax(lambda),xstyle=1
   oplot,lambda,cont,color=fsc_color('yellow')
   contdiv = spec/cont
   
   ;fit for age and velocity dispersion
   ;fix metallicity to -0.4 dex
    sps = sps_read_spec('/scr2/nichal/workspace2/fsps-3.0/SSP/SSP_Padova_MILES_Kroupa_zmet16.out.spec')
    spsage = sps.agegyr
    a = [3.,10.,redshift] ;age,vdisp,z
    fit = curvefit(lambda, contdiv,weight,a,function_name='quickspec')

   plot,lambda,contdiv
   oplot,lambda,fit,color=fsc_color('red')
 
   ;calculate signal to noise
   dev = abs((contdiv-fit)/fit)
   avgdev = mean(dev)
   w = where(dev lt 3.0*avgdev, c)
   if c gt 0 then sn_perpix = 1.0/mean(dev[w])

   sn = sn_perpix*sqrt(1./cdelt)
   print,'S/N per pix,angstrom from fit is ',sn_perpix,sn

   ;signal to noise from continuum
   dev = abs(contdiv-1.)
   avgdev = mean(dev)
   w = where(dev lt 3.0*avgdev, c)
   if c gt 0 then sn_perpix = 1.0/mean(dev[w])
   sn = sn_perpix*sqrt(1./cdelt)
   print,'S/N per pix,angstrom from continuum is ',sn_perpix,sn

   stop
end
