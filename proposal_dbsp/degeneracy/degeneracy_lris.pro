pro degeneracy_lris,sn=sn,redo=redo,fine=fine
common degen, spec_Arr, grid_feh, grid_age
	if keyword_set(redo) eq 1 and keyword_set(fine) eq 0 then begin
		;make spectra with these grids
		grid_feh = findgen(121)/100.-1.2 ;range from [-1.2,0.0] dex
		grid_age = findgen(126)/10.+0.5 ;range from [0.5,13] Gyr
		spec_arr = make_spec_lris(grid_feh,grid_age)
		mwrfits,spec_arr,'grid_spec_lris.fits',/silent,/create
	endif 
	if keyword_set(redo) eq 1 and keyword_set(fine) eq 1 then begin
                ;make spectra with these grids
                grid_feh = findgen(601)/500.-1.2 ;range from [-1.2,0.0] dex
                grid_age = findgen(12)+0.5 ;range from [0.5,11.5] Gyr
                spec_arr = make_spec_lris(grid_feh,grid_age)
                mwrfits,spec_arr,'grid_spec_lris_fine.fits',/silent,/create		
        endif
	if keyword_set(fine) then spec_arr = mrdfits('grid_spec_lris_fine.fits',1) else $
		 spec_arr = mrdfits('grid_spec_lris.fits',1)
	fehpos = uniq(spec_arr.feh,sort(Spec_arr.feh))
	grid_feh = spec_Arr(fehpos).feh
	agepos = uniq(spec_arr.age,sort(Spec_arr.age))
	grid_age = spec_Arr(agepos).age

	;find chisquare at each reference point
	if keyword_set(fine) eq 0 then begin
		ref_ia = [20,30,40,50,60,70] ;2,5,3.5,4.5,5.5,6.5,7.5 Gyr
		ref_iz = [20,30,40,50,60,70,80,90,100,110] ;-1,-0.9,-0.8,...-0.1
        endif else begin
                ref_ia = [2,3,4,5,6,7] ;2,5,3.5,4.5,5.5,6.5,7.5 Gyr
                ref_iz = [20,30,40,50,60,70,80,90,100,110]*5 ;-1,-0.9,-0.8,...-0.1
	endelse
	for i=0,n_Elements(ref_ia)-1 do for j=0,n_elements(ref_iz)-1 do begin
		ageref = grid_age(ref_ia(i))
		fehref = grid_feh(ref_iz(j))
		uncert = make_chisqarr_lris(ref_ia(i),ref_iz(j),sn)
	endfor
	;stop
end
