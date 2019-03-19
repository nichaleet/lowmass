function make_spec_lris,grid_feh,grid_age
common sps_spec, sps, spsz, spsage
common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize,rest
	;setting
        ;for 600/4000 grating for blue side
	redshift = 0.02
	vdisp = 90.  ;FWHM (approximately for 10^9 Msun in Dutton11a)
        lambdamin = 3300.
        dlam = 4./2.35 ;angstrom FWHM/2.35
	npix = 4096
        dispersion = 0.63 ;angstrom per pixel
        
	datalam = lambdamin+findgen(npix)*dispersion 
	dlam = dlam+fltarr(npix)
 
	wonfit = 1+bytarr(npix)

	xmp = datalam(wonfit)
	normalize = 0
	rest = 0
	
	;loop over grid_feh and grid_age
	nfeh = n_elements(grid_feh)
	nage = n_elements(grid_age)
	for iz=0,nfeh-1 do for ia=0,nage-1 do begin
		spsspec = get_sps_obs(xmp,[grid_feh(iz),grid_age(ia),vdisp/2.35,redshift])
		str = {lambda:xmp,spec:spsspec,feh:grid_feh(iz),age:grid_age(ia),vdisp:vdisp,redshift:redshift}	
		if iz eq 0 and ia eq 0 then strall=replicate(str,nfeh,nage)
		strall[iz,ia] = str
		if iz mod 10 eq 0 and ia mod 10 eq 0 then print, iz, ia
	endfor
	return,strall
end
