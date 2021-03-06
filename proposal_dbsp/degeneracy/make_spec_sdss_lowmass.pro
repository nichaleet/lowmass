function make_spec_sdss_lowmass,grid_feh,grid_age
common sps_spec, sps, spsz, spsage
common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize,rest
	;setting
	redshift = 0.05
	vdisp = 90./2.35

	;get the obs setting	
	sdss = mrdfits('/scr2/nichal/workspace2/sps_fit/data/gallazzi_allmass2/sps_fit.fits.gz',1)
	science = sdss[0]
	reallambda = science.lambda 
	dlam = science.dlam
 
	znow = redshift
	won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and reallambda/(1.+znow) gt 3500. and reallambda/(1.+znow) lt 5900., con)
	wonfit = won

	xmp = reallambda[won]
	ymp = science.contdiv[won]
	dymp = (science.contdivivar[won])^(-0.5)

	dataivar = science.contdivivar
	datalam = reallambda
	contmask = science.contmask
	normalize = 0
	rest = 0
	
	;loop over grid_feh and grid_age
	nfeh = n_elements(grid_feh)
	nage = n_elements(grid_age)
	for iz=0,nfeh-1 do for ia=0,nage-1 do begin
		spsspec = get_sps_obs_sdss(xmp,[grid_feh(iz),grid_age(ia),vdisp,redshift])
		str = {lambda:xmp,spec:spsspec,feh:grid_feh(iz),age:grid_age(ia),vdisp:vdisp,redshift:redshift}		
		if iz eq 0 and ia eq 0 then strall=replicate(str,nfeh,nage)
		strall[iz,ia] = str
		if iz mod 10 eq 0 and ia mod 10 eq 0 then print, iz, ia
	endfor
	return,strall
end
