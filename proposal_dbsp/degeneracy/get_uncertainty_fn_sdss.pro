pro get_uncertainty_fn_sdss
;to get theoretical uncertainties as a function of Fe/H and SN
firstit = 1
set_plot,'x'
device,decomposed=0
loadct,23
!p.multi=[0,1,2]
;loop over signal to noise
for iu=10,52,1 do begin
	dir='sdsslowmassoutputsn'+strtrim(string(fix(iu)),2)
	;dir='outputsn'+strtrim(string(fix(iu)),2)
	print, 'doing '+dir
	files = file_search(dir+'/chisq*.fits')
	;get values from each file
	for i=0,n_Elements(files)-1 do begin
		;get age and metallicities from file name
		words = strsplit(files(i),'_',/extract)
		age = float(strmid(words[1],3))
		feh = float(strmid(words[2],3,strpos(words[2],'.fits')))
		;read the file
		chisqarr = mrdfits(files(i),0,header)
		;you should check below if it matches those in degeneracy.pro
		grid_Feh = rebin(indgen(601)/500.-1.2,601,12)
		grid_age = transpose(rebin(alog10(findgen(12)+0.5)+9.,12,601)) ;range from [0.5,13] Gyr		
		minchisq = min(chisqarr)
		goodchisq = where(chisqarr lt minchisq+1.,ngoodchisq)
		if ngoodchisq lt 2 then stop
		agerange = minmax(grid_age(goodchisq))
		fehrange = minmax(grid_feh(goodchisq))
		str = {age:age,feh:feh,agerange:agerange,fehrange:fehrange,bestage:grid_age(minchisq),bestfeh:grid_feh(minchisq)}
		if i eq 0 then fixsnstr = str else fixsnstr=[fixsnstr,str]	
	endfor
	feharr = fixsnstr(uniq(fixsnstr.feh,sort(fixsnstr.feh))).feh
	nfeh = n_Elements(feharr)
	sigmaage = fltarr(nfeh) 
	sigmafeh_upper = fltarr(nfeh)
	sigmafeh_lower = fltarr(nfeh)
	;run loop over each fixed metallicity
	for i=0,nfeh-1 do begin
		;choose all str with the metallicities, multiple ages)
		selfeh = where(fixsnstr.feh eq feharr(i),cselfeh)
		delage = 0.5*(fixsnstr(selfeh).agerange[1]-fixsnstr(selfeh).agerange[0])
		delfeh_upper = abs(fixsnstr(selfeh).fehrange[1]-fixsnstr(selfeh).feh)
		delfeh_lower = abs(fixsnstr(selfeh).fehrange[0]-fixsnstr(selfeh).feh)
		sigmaage(i) = median(delage) ;in dex
		sigmafeh_upper(i) = median(delfeh_upper) ;in dex
		sigmafeh_lower(i) = median(delfeh_lower) ;in dex
	endfor
	if firstit eq 1 then begin
		plot,feharr,sigmaage,xtitle='[Fe/H]',ytitle='sigma log(age)',color=fsc_color('white'),title='SDSS SN=[20,40]'
		plot,feharr,sigmafeh_upper,xtitle='[Fe/H]',ytitle='sigma FEH',yrange=[0,0.2],color=fsc_color('white')
	endif else begin
		oplot,feharr,sigmafeh_upper,color=(iu-19)*7,linestyle=2
		oplot,feharr,sigmafeh_lower,color=(iu-19)*7
	endelse
	outstr = {sn:iu,feh:feharr,sigmaage:sigmaage,sigmafeh_upper:sigmafeh_upper,sigmafeh_lower:sigmafeh_lower}
	if firstit eq 1 then finstr = outstr else finstr=[finstr,outstr]
	firstit = 0
endfor

mwrfits,finstr,'sdss_lowmass_uncertainty_fn_feh_sn.fits',/silent,/create

end
