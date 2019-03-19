pro get_uncertainty_fn_lris

;to get theoretical uncertainties as a function of Fe/H and SN
firstit = 1
set_plot,'x'
device,decomposed=0
loadct,23
!p.multi=[0,1,2]

;loop over signal to noise
for iu=20,52 do begin
	dir='lrisoutputsn'+strtrim(string(fix(iu)),2)
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
;		if i le 39 then begin
;			grid_Feh = rebin(findgen(121)/100.-1.2,121,126) 
;			grid_age = transpose(rebin(alog10(findgen(126)/10.+0.5)+9.,126,121)) ;range from [0.5,13] Gyr	
;		endif else begin
			grid_feh = rebin(findgen(601)/500.-1.2,601,12)
			grid_age = transpose(rebin(alog10(findgen(12)+0.5)+9.,12,601))
;		endelse	
		minchisq = min(chisqarr)
		goodchisq = where(chisqarr lt minchisq+1.,ngoodchisq)
		if ngoodchisq lt 2 then stop
		agerange = minmax(grid_age(goodchisq))
		fehrange = minmax(grid_feh(goodchisq))
		str = {age:age,feh:feh,agerange:agerange,fehrange:fehrange,bestfeh:grid_feh(minchisq),bestage:grid_age(minchisq)}
		if i eq 0 then fixsnstr = str else fixsnstr=[fixsnstr,str]	
	endfor
	feharr = fixsnstr(uniq(fixsnstr.feh,sort(fixsnstr.feh))).feh
	nfeh = n_Elements(feharr)
	sigmaage = fltarr(nfeh) 
	sigmafeh = fltarr(nfeh)
	;run loop over each fixed metallicity
	for i=0,nfeh-1 do begin
		;choose all str with the metallicities, multiple ages)
		selfeh = where(fixsnstr.feh eq feharr(i),cselfeh)
		delage = 0.5*(fixsnstr(selfeh).agerange[1]-fixsnstr(selfeh).agerange[0])
		delfeh = 0.5*(fixsnstr(selfeh).fehrange[1]-fixsnstr(selfeh).fehrange[0])
		sigmaage(i) = median(delage) ;in dex
		sigmafeh(i) = median(delfeh) ;in dex
	endfor
	if firstit eq 1 then begin
                plot,feharr,sigmaage,xtitle='[Fe/H]',ytitle='sigma log(age)',color=fsc_color('white'),title='LRIS SN=[20,52]'
                plot,feharr,sigmafeh,xtitle='[Fe/H]',ytitle='sigma FEH',yrange=[0,0.5],color=fsc_color('white')
        endif else begin
                oplot,feharr,sigmaage,color=(iu-19)*10
                oplot,feharr,sigmafeh,color=(iu-19)*10
        endelse

	outstr = {sn:iu,feh:feharr,sigmaage:sigmaage,sigmafeh:sigmafeh}
        if firstit eq 1 then finstr = outstr else finstr=[finstr,outstr]
        firstit = 0
	;stop
endfor

mwrfits,finstr,'lris_uncertainty_fn_feh_sn.fits',/silent,/create
end
