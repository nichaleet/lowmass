pro testslope_cursdss
   sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/lowmass/sps_fit.fits.gz',1)
   cat = mrdfits('/scr2/nichal/workspace2/lowmass/ana/lowmass_nichaleet_match.fits',1)
   good = where(sci.good eq 1 and sci.feherr gt 0. and sci.goodfit eq 0 and cat.logmass_min gt 0. and cat.logmass_max gt 0.,cgood)
   mass = sci(good).logmstar
   sn = cat(good).snmedian
   feherr = sci(good).feherr

   sci2 = mrdfits('/scr2/nichal/workspace2/sps_fit/data/lowmass_lowsn/sps_fit.fits.gz',1)
   cat2 = mrdfits('/scr2/nichal/workspace2/lowmass/ana/lowmass_lowsn_nichaleet_match.fits',1)
   good2 = where(sci2.good eq 1 and sci2.feherr gt 0. and sci2.goodfit eq 0 and cat2.logmass_min gt 0. and cat2.logmass_max gt 0.,cgood2)
   mass2 = sci2(good2).logmstar
   sn2 = cat2(good2).snmedian
   feherr2 = sci2(good2).feherr

   mass = [mass,mass2]
   feherr = [feherr,feherr2]
   ngals = n_elements(mass)
   sn = [sn,sn2]

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
   intrinsic_scat = 0.1
   brd_slope = 0.95/1.5  ;brd = bridge
   brd_const = -0.45 ;constant at 10^9 Msun
   brd_rms = 0.1
   frac_enk = 0.3 ; fractions of gals that follow dwarf galaxies

   marr = findgen(16)/10.+8.
   enk_line = enk_slope*(marr-6)+enk_const
   nl_line = nl_slope*(marr-10)+nl_const
   brd_line = brd_slope*(marr-9)+brd_const

   mass_Brd = mass
   ngals_enk = round(frac_enk*ngals)
   ngals_nl = ngals-ngals_enk
   ;randomly assign which one is enk which one is nl
   ord = sort(randomu(seed,ngals))
   mass_enk = mass(ord[0:ngals_enk-1])
   mass_nl  = mass(ord[ngals_enk:ngals-1])
   feherr_enk = feherr(ord[0:ngals_enk-1])
   feherr_nl = feherr(ord[ngals_enk:ngals-1])
   sn_enk = randomn(seed,ngals_enk)*8+14.
   sn_nl = randomn(seed,ngals_nl)*8+14
   sn_brd = randomn(seed,ngals)*8+14.

   niter = 2
   iterplot = 1
   outstr = {D_ks:fltarr(niter),pvalue:fltarr(niter)}

   for iiter=0,niter-1 do begin
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

        
           deltaobsfeh_upper_enk = sqrt(deltaobsfeh_upper_enk^2+feherr_enk^2)
           deltaobsfeh_lower_enk = sqrt(deltaobsfeh_lower_enk^2+feherr_enk^2)
           deltaobsfeh_upper_nl = sqrt(deltaobsfeh_upper_nl^2+feherr_nl^2)
           deltaobsfeh_lower_nl = sqrt(deltaobsfeh_lower_nl^2+feherr_nl^2)
           set_plot,'ps'
           psname='cursdss_example_mzr_dbsp_case1.eps'
           device,filename=psname,xsize=15,ysize=10,/encapsulated,/color
           !p.font = 0
           !p.multi=[0,1,1]
	   ymin = -1.2
	   plot,[mass,mass],[obs_feh_nl,obs_feh_enk],xtitle='log(mass)',ytitle='[Fe/H]',psym=1,/nodata,$
                   yrange=[ymin,0.25],ystyle=1,title='Simulation with current SDSS S/Ns'
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

           vsym,4,rot=45,/fill
           oplot,mass_nl,obs_feh_nl,psym=8,color=fsc_color('red')
           oplot,mass_enk,obs_feh_enk,psym=8,color=fsc_color('navy')
	   xyouts,8.1,0.15,'Case I: bimodal MZR'
           al_legend,['Kirby et al.,2013','Leethochawalit et al., in prep'],psym=[14,14],color=fsc_color(['navy','red']),position=[8.0,0.15],box=0
           device,/close

           psname='cursdss_example_mzr_dbsp_case2.eps'
           device,filename=psname,xsize=15,ysize=10,/encapsulated,/color
           !p.font = 0
           !p.multi=[0,1,1]
	   ymin = -1.2
	   plot,[mass,mass],[obs_feh_nl,obs_feh_enk],xtitle='log(mass)',ytitle='[Fe/H]',psym=1,/nodata,$
                   yrange=[ymin,0.25],ystyle=1,title='Simulation with current SDSS S/Ns'
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

           cgerrplot,mass_brd,obs_feh_brd-deltaobsfeh_lower_brd,obs_feh_brd+deltaobsfeh_upper_brd

           vsym,4,rot=45,/fill
           oplot,mass_brd,obs_Feh_brd,psym=cgsymcat(16),color=fsc_color('darkorchid')
	   xyouts,8.1,0.15,'Case II: continuous MZR'
         ;  al_legend,['Case I: bimodal MZR','Case II: continuous MZR'],psym=[14,16],color=fsc_color(['red','darkorchid']),/top,/left,box=0
          ; oplot,[8.05],[0.115],psym=cgsymcat(14),color=fsc_color('navy')
           device,/close
	endif
	ks_test,mass,obs_Feh_enk,mass,obs_feh_nl,D_ks,Neff_ks,pvalue,/silent,/noplot
        outstr.d_ks(iiter) = d_ks
        outstr.pvalue(iiter) = pvalue
   endfor
 
   if niter gt 100 then begin
      psname='cursdss_ks_param_dbsp.eps'
      device, filename = psname,xsize = 15,ysize = 20, $
         xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      !p.multi=[0,1,2]
      !p.font = 0
      !p.charsize=0.5
      plothist,outstr.d_ks,xtitle='KS Distance'
      plothist,outstr.pvalue,xtitle='pvalue',bin=0.05
      device,/close
   endif 
end
