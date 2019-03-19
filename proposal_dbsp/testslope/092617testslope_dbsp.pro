pro testslope_dbsp
;double spec
enk_slope = 0.3   ;dex
enk_slope_err = 0.02 ;dex
enk_rms = 0.17    ;dex
enk_const = -1.69 ;dex this is constant at 10^6 Msun
enk_const_err = 0.04 ;dex
nl_slope = 0.15   ;dex
nl_slope_err = 0.03 ;dex
nl_const = -0.07 ;dex this is constant at 10^10 msun (y=m(x-10)+c)
nl_const_err = 0.01 ;dex

mag = 18.5
n_nights = 1
hrs_per_night = 5
;600 line grating(blue) center at 5000 will cover from 3500 to 6500, R=1806

niter = 1000
iterplot = 500.

itime_arr = [4]
outstr = {itime:5.,sn:0.,ngals:0,D_ks:fltarr(niter),pvalue:fltarr(niter)}
outstr = replicate(outstr,n_elements(itime_arr))
for i=0,n_elements(itime_arr)-1 do begin
   itime = itime_arr(i)
   outstr[i].itime = itime
   ngals = hrs_per_night*60/itime*n_nights
   sn = fltarr(ngals)+calsn_dbsp(mag,itime*60,5000.,600) 
   outstr[i].ngals=ngals
   print, 'int time',itime, ' ngals=',ngals, ' sn=',sn[0]
   for iiter=0,niter-1 do begin
   	mass = randomu(seed,ngals)+8.5
   	feh_enk = enk_slope*(mass-6)+enk_const;+0.57 ;add 0.57 so that the FeH at 10^9 Msun is the same for both enk and nl
   	feh_nl = nl_slope*(mass-10)+nl_const

   	get_additional_uncertainties_lowmass,feh_enk,sn,deltaage_enk,deltafeh_upper_enk,deltafeh_lower_enk,/sdss
   	get_additional_uncertainties_lowmass,feh_nl,sn,deltaage_nl,deltafeh_upper_nl,deltafeh_lower_nl,/sdss

	obs_feh_enk_upper = feh_enk+abs(randomn(seed,ngals)*deltafeh_upper_enk)
	obs_feh_enk_lower = feh_enk-abs(randomn(seed,ngals)*deltafeh_lower_enk)
	low_where = where(randomu(seed,ngals) lt 0.5, complement=high_where)
	obs_feh_enk = fltarr(ngals)
        obs_feh_enk(low_where) = obs_feh_enk_lower(low_where)
	obs_feh_enk(high_where) = obs_feh_enk_upper(high_where)

	obs_feh_nl_upper = feh_nl+abs(randomn(seed,ngals)*deltafeh_upper_nl)
        obs_feh_nl_lower = feh_nl-abs(randomn(seed,ngals)*deltafeh_lower_nl)
        low_where = where(randomu(seed,ngals) lt 0.5, complement=high_where)
        obs_feh_nl = fltarr(ngals)
        obs_feh_nl(low_where) = obs_feh_nl_lower(low_where)
        obs_feh_nl(high_where) = obs_feh_nl_upper(high_where)

	get_additional_uncertainties_lowmass,obs_feh_enk,sn,deltaobsage_enk,deltaobsfeh_upper_enk,deltaobsfeh_lower_enk,/sdss
        get_additional_uncertainties_lowmass,obs_feh_nl,sn,deltaobsage_nl,deltaobsfeh_upper_nl,deltaobsfeh_lower_nl,/sdss
        	
        if iiter eq iterplot then begin
           set_plot,'ps'
           psname='int_time'+sigfig(itime,2)+'mins_'+sigfig(n_nights,2)+'nights_example_mzr_dbsp.eps'
           device,filename=psname,xsize=15,ysize=10,/encapsulated,/color
	   !p.font = 0
           !p.multi=[0,1,1]
           plot,[mass,mass],[obs_feh_nl,obs_feh_enk],xtitle='logmass',ytitle='[Fe/H]',psym=1/nodata
           cgerrplot,mass,obs_feh_nl-deltaobsfeh_lower_nl,obs_feh_nl+deltaobsfeh_upper_nl
           oploterror,mass,obs_feh_enk,deltaobsfeh_enk,psym=1,errocolor=fsc_color('navy')
	   vsym,4,rot=45,/fill
           oplot,mass,obs_feh_nl,psym=8,color=fsc_color('red')
           oplot,mass,obs_feh_enk,psym=8,color=fsc_color('navy')
           device,/close
	print, 'mean delta feh',mean(deltafeh_enk),mean(deltafeh_nl)
        endif
           
	ks_test,mass,obs_Feh_enk,mass,obs_feh_nl,D_ks,Neff_ks,pvalue,/silent,/noplot	
	outstr[i].d_ks(iiter) = d_ks
        outstr[i].pvalue(iiter) = pvalue
   endfor
   psname='int_time'+sigfig(itime,2)+'mins_'+sigfig(n_nights,2)+'nights_ks_param_dbsp.eps'
   device, filename = psname,xsize = 15,ysize = 20, $
       xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
   !p.multi=[0,1,2]
   !p.font = 0
   plothist,outstr[i].d_ks,xtitle='KS Distance'
   plothist,outstr[i].pvalue,xtitle='pvalue'
   device,/close

endfor

save,outstr,filename='testslope_dbsp.sav'
 
end
