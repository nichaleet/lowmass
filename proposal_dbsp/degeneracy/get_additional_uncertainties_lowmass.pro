pro get_additional_uncertainties_lowmass,feh,sn,deltaage,deltafeh_upper,deltafeh_lower,lris=lris,sdss=sdss
        ;get the uncertainties files
	
        if keyword_set(lris) then delstr = mrdfits('/scr2/nichal/workspace2/lowmass/degeneracy/uncertainty_fn_feh_sn_lris.fits',1,/silent)
	if keyword_set(sdss) then delstr = mrdfits('/scr2/nichal/workspace2/lowmass/degeneracy/sdss_lowmass_uncertainty_fn_feh_sn.fits',1,/silent)
        nobjs = n_Elements(feh)
        deltaage = fltarr(nobjs);in dex
        deltafeh_upper = fltarr(nobjs)
        deltafeh_lower = fltarr(nobjs)
        for i=0,nobjs-1 do begin
           boobee = min(abs(delstr.sn-sn(i)),loc)
           deltaage(i) = interpol(delstr(loc).sigmaage,delstr(loc).feh,feh(i))
           deltafeh_upper(i) = interpol(delstr(loc).sigmafeh_upper,delstr(loc).feh,feh(i))
           deltafeh_lower(i) = interpol(delstr(loc).sigmafeh_lower,delstr(loc).feh,feh(i))
        endfor
;stop
end
~        
