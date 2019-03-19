pro degeneracy_sdss_lowmass,sn=sn,redo=redo
;e.g. for i=20,40,1 do degeneracy_sdss_lowmass,sn=i,/redo
common degen, spec_Arr, grid_feh, grid_age
	if n_Elements(redo) eq 1 then begin
		;make spectra with these grids
		grid_feh = findgen(601)/500.-1.2 ;range from [-1.2,0.] dex 
		grid_age = findgen(12)+0.5 ;range from [0.5,13] Gyr
		spec_arr = make_spec_sdss(grid_feh,grid_age)
		mwrfits,spec_arr,'grid_spec_sdss_lowmass.fits',/silent,/create	
	endif else spec_arr = mrdfits('grid_spec_sdss_lowmass.fits',1)
	fehpos = uniq(spec_arr.feh,sort(Spec_arr.feh))
	grid_feh = spec_Arr(fehpos).feh
	agepos = uniq(spec_arr.age,sort(Spec_arr.age))
	grid_age = spec_Arr(agepos).age

	;find chisquare at each reference point
	ref_ia = [2,3,4,5,6,7] ;2,5,3.5,4.5,5.5,6.5,7.5 Gyr
	ref_iz = [20,30,40,50,60,70,80,90,100,110]*5 ; -0.44,-0.36,-0.28,-0.2,-0.12,-0.04,0.04,0.12
	for i=0,n_Elements(ref_ia)-1 do for j=0,n_elements(ref_iz)-1 do begin
		ageref = grid_age(ref_ia(i))
		fehref = grid_feh(ref_iz(j))
		uncert = make_chisqarr_sdss_lowmass(ref_ia(i),ref_iz(j),sn)
	endfor
end
