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
nl_rms = 0.12
brd_slope = 0.95/1.5  ;brd = bridge
brd_const = -0.45 ;constant at 10^9 Msun
brd_rms = 0.1 
frac_enk = 0.3 ; fractions of gals that follow dwarf galaxies

intrinsic_scat = 0.1 ;dex

itime_arr = [15,20,30,30]
n_nights = [3,3,3,4]
hrs_per_night = 6.5
;600 line grating(blue) center at 5000 will cover from 3500 to 6500, R=1806


marr = findgen(16)/10.+8.
enk_line = enk_slope*(marr-6)+enk_const
nl_line = nl_slope*(marr-10)+nl_const
brd_line = brd_slope*(marr-9)+brd_const

niter = 2
iterplot = 1.

outstr = {itime:5.,sn:0.,ngals:0,D_ks:fltarr(niter),pvalue:fltarr(niter)}
outstr = replicate(outstr,n_elements(itime_arr))
for i=0,n_elements(itime_arr)-1 do begin
   itime = itime_arr(i)
   n_night = n_nights[i]
   outstr[i].itime = itime
   ngals = round(hrs_per_night*60/itime*n_night)
   ngals_enk = round(ngals*frac_enk)
   ngals_nl = ngals-ngals_enk

   mag_enk = fltarr(ngals_enk)+randomu(seed,ngals_enk)*1.5+17.5 ;random 17.5 to 19 mag
   sn_enk = fltarr(ngals_enk)+calsn_dbsp(mag_enk,itime*60,5000.,600) 
   mag_nl = fltarr(ngals_nl)+randomu(seed,ngals_nl)*1.5+17.5 ;random 17.5 to 19 mag
   sn_nl = fltarr(ngals_nl)+calsn_dbsp(mag_nl,itime*60,5000.,600)
   mag_brd = fltarr(ngals)+randomu(seed,ngals)*1.5+17.5    
   sn_brd = fltarr(ngals)+calsn_dbsp(mag_brd,itime*60,5000.,600)

   outstr[i].ngals=ngals
   print, 'int time',itime, 'min ngals(bridge)=',ngals, ' sn(mean,min,max)=',mean(sn_brd),minmax(sn_brd)
   print, 'int time',itime, 'min ngals(enk)=',ngals_enk, ' sn(mean,min,max)=',mean(sn_enk),minmax(sn_enk)
   print, 'int time',itime, 'min ngals(nl)=',ngals_nl, ' sn(mean,min,max)=',mean(sn_nl),minmax(sn_nl)

   for iiter=0,niter-1 do begin
   	mass_enk = randomu(seed,ngals_enk)*1.5+8.
	mass_nl  = randomu(seed,ngals_nl)*1.5+8.
        mass_brd = randomu(seed,ngals)*1.5+8.
   	feh_enk = enk_slope*(mass_enk-6)+enk_const+randomu(seed,ngals_enk)*intrinsic_scat
   	feh_nl = nl_slope*(mass_nl-10)+nl_const+randomu(seed,ngals_nl)*intrinsic_scat
	feh_brd = brd_slope*(mass_brd-9)+brd_const+randomu(seed,ngals)*intrinsic_scat

   	get_additional_uncertainties_lowmass,feh_enk,sn_enk,deltaage_enk,deltafeh_upper_enk,deltafeh_lower_enk,/sdss
   	get_additional_uncertainties_lowmass,feh_nl,sn_nl,deltaage_nl,deltafeh_upper_nl,deltafeh_lower_nl,/sdss
        get_additional_uncertainties_lowmass,feh_brd,sn_brd,deltaage_brd,deltafeh_upper_brd,deltafeh_lower_brd,/sdss

        obs_feh_enk_upper = feh_enk+abs(randomn(seed,ngals_enk)*deltafeh_upper_enk)
        obs_feh_enk_lower = feh_enk-abs(randomn(seed,ngals_enk)*deltafeh_lower_enk)
        low_where = where(randomu(seed,ngals_enk) lt 0.5, complement=high_where)
        obs_feh_enk = fltarr(ngals_enk)
        obs_feh_enk(low_where) = obs_feh_enk_lower(low_where)
        obs_feh_enk(high_where) = obs_feh_enk_upper(high_where)

        obs_feh_nl_upper = feh_nl+abs(randomn(seed,ngals_nl)*deltafeh_upper_nl)
        obs_feh_nl_lower = feh_nl-abs(randomn(seed,ngals_nl)*deltafeh_lower_nl)
        low_where = where(randomu(seed,ngals_nl) lt 0.5, complement=high_where)
        obs_feh_nl = fltarr(ngals_nl)
        obs_feh_nl(low_where) = obs_feh_nl_lower(low_where)
        obs_feh_nl(high_where) = obs_feh_nl_upper(high_where)

        obs_feh_brd_upper = feh_brd+abs(randomn(seed,ngals)*deltafeh_upper_brd)
        obs_feh_brd_lower = feh_brd-abs(randomn(seed,ngals)*deltafeh_lower_brd)
        low_where = where(randomu(seed,ngals) lt 0.5, complement=high_where)
        obs_feh_brd = fltarr(ngals)
        obs_feh_brd(low_where) = obs_feh_brd_lower(low_where)
        obs_feh_brd(high_where) = obs_feh_brd_upper(high_where)
        	
        if iiter eq iterplot then begin

	   get_additional_uncertainties_lowmass,obs_feh_enk,sn_enk,deltaobsage_enk,deltaobsfeh_upper_enk,deltaobsfeh_lower_enk,/sdss
           get_additional_uncertainties_lowmass,obs_feh_nl,sn_nl,deltaobsage_nl,deltaobsfeh_upper_nl,deltaobsfeh_lower_nl,/sdss
           get_additional_uncertainties_lowmass,obs_feh_brd,sn_brd,deltaobsage_brd,deltaobsfeh_upper_brd,deltaobsfeh_lower_brd,/sdss

           set_plot,'ps'
           psname='int_time'+sigfig(itime,2)+'mins_'+sigfig(n_night,2)+'nights_example_mzr_dbsp.eps'
           device,filename=psname,xsize=15,ysize=10,/encapsulated,/color
	   !p.font = 0
           !p.multi=[0,1,1]
	   ymin = -1.2
           plot,[mass_nl,mass_enk],[obs_feh_nl,obs_feh_enk],xtitle='log(mass)',ytitle='[Fe/H]',psym=1,/nodata,$
                     yrange=[ymin,0.25],ystyle=1,xrange=[8,9.5],xstyle=1,title='Simulation with DBSP'
           polyfill,[marr,reverse(marr)],[(enk_line+enk_rms)>(ymin),(reverse(enk_line)-enk_rms)>(ymin)],color=fsc_color('blu2')
           oplot,marr,enk_line,linestyle=2,color=fsc_color('blu5')
           polyfill,[marr,reverse(marr)],[(nl_line+nl_rms)>(ymin),(reverse(nl_line)-nl_rms)>(ymin)],color=fsc_color('org2')
           oplot,marr,nl_line,linestyle=2,color=fsc_color('org5')
           polyfill,[marr,reverse(marr)],[(brd_line+brd_rms)>(ymin),(reverse(brd_line)-brd_rms)>(ymin)],color=fsc_color('thistle'),/line_fill,orientation=45
           polyfill,[marr,reverse(marr)],[(brd_line+brd_rms)>(ymin),(reverse(brd_line)-brd_rms)>(ymin)],color=fsc_color('thistle'),/line_fill,orientation=135
           oplot,marr,brd_line,linestyle=2,color=fsc_color('purple')
            
           axis,yaxis=0,yrange=[ymin,0.25],ystyle=1,ytickformat='(A1)'
           axis,yaxis=1,yrange=[ymin,0.25],ystyle=1,ytickformat='(A1)'
           axis,xaxis=0,xrange=[8,9.5],xtickformat='(A1)'

           cgerrplot,mass_nl,obs_feh_nl-deltaobsfeh_lower_nl,obs_feh_nl+deltaobsfeh_upper_nl
	   cgerrplot,mass_enk,obs_feh_enk-deltaobsfeh_lower_enk,obs_feh_enk+deltaobsfeh_upper_enk
	   cgerrplot,mass_brd,obs_feh_brd-deltaobsfeh_lower_brd,obs_feh_brd+deltaobsfeh_upper_brd

	   vsym,4,rot=45,/fill
           oplot,mass_nl,obs_feh_nl,psym=8,color=fsc_color('red')
           oplot,mass_enk,obs_feh_enk,psym=8,color=fsc_color('navy')
	   oplot,mass_brd,obs_Feh_brd,psym=cgsymcat(16),color=fsc_color('darkorchid')
	   al_legend,['Case I: bimodal MZR','Case II: continuous MZR'],psym=[14,16],color=fsc_color(['red','darkorchid']),/top,/left,box=0
	   oplot,[8.05],[0.115],psym=cgsymcat(14),color=fsc_color('navy')
           device,/close
	print, 'mean delta feh upper',mean(deltafeh_upper_enk),mean(deltafeh_upper_nl)
        print, 'mean delta feh lower',mean(deltafeh_lower_enk),mean(deltafeh_lower_nl)

        endif
           
	ks_test,[mass_enk,mass_nl],[obs_Feh_enk,obs_feh_nl],mass_brd,obs_feh_brd,D_ks,Neff_ks,pvalue,/silent,/noplot	
	outstr[i].d_ks(iiter) = d_ks
        outstr[i].pvalue(iiter) = pvalue
   endfor

   if niter gt 100 then begin
      psname='int_time'+sigfig(itime,2)+'mins_'+sigfig(n_night,2)+'nights_ks_param_dbsp.eps'
      device, filename = psname,xsize = 15,ysize = 20, $
          xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.multi=[0,1,2]
      !p.font = 0
      !p.charsize=0.5
      plothist,outstr[i].d_ks,xtitle='KS Distance'
      plothist,outstr[i].pvalue,xtitle='pvalue',bin=0.05
   endif
   device,/close

endfor

save,outstr,filename='testslope_dbsp.sav'
 
end
