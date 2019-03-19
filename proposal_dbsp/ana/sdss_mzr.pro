pro sdss_mzr

   enk_slope = 0.3   ;dex
   enk_slope_err = 0.02 ;dex
   enk_rms = 0.17    ;dex
   enk_const = -1.69 ;dex this is constant at 10^6 Msun
   enk_const_err = 0.04 ;dex
   nl_rms = 0.13
   nl_slope = 0.15   ;dex
   nl_slope_err = 0.03 ;dex
   nl_const = -0.07 ;dex this is constant at 10^10 msun (y=m(x-10)+c)
   nl_const_err = 0.01 ;dex

   ;lowmass stuff
   sci = mrdfits('/scr2/nichal/workspace2/sps_fit/data/lowmass/sps_fit.fits.gz',1)
   cat = mrdfits('/scr2/nichal/workspace2/lowmass/ana/lowmass_nichaleet_match.fits',1)
   restore,'lowmass_Ha_ew.sav'
   good = where(sci.good eq 1 and sci.feherr gt 0. and sci.goodfit eq 0 and cat.logmass_min gt 0. and cat.logmass_max gt 0.,cgood)
   print,'final sample size=',cgood
   feh = sci(good).feh
   mass = cat(good).logmass
   badmassmin = where(cat(good).logmass_min lt 0.,cbadmassmin)
   badmassmax = where(cat(good).logmass_max lt 0.,cbadmassmax)
   massmin = cat(good).logmass_min
   massmax = cat(good).logmass_max
   massmin(badmassmin) = mass(badmassmin)
   massmax(badmassmax) = mass(badmassmax)
   feherr = sci(good).feherr
   sn = sci(good).sn
   get_additional_uncertainties_lowmass,feh,sn,deltaage,deltafeh_upper,deltafeh_lower,/sdss
   feherr_upper = sqrt(feherr^2+deltafeh_upper^2)
   feherr_lower = sqrt(feherr^2+deltafeh_lower^2)

   ;lowsn stuff
   sci2 = mrdfits('/scr2/nichal/workspace2/sps_fit/data/lowmass_lowsn/sps_fit.fits.gz',1)
   cat2 = mrdfits('/scr2/nichal/workspace2/lowmass/ana/lowmass_lowsn_nichaleet_match.fits',1)
   restore,'lowmass_lowsn_Ha_ew.sav'
   good2 = where(sci2.good eq 1 and sci2.feherr gt 0. and sci2.goodfit eq 0 and cat2.logmass_min gt 0. and cat2.logmass_max gt 0.,cgood2)
   print,'final sample size2=',cgood2
   feh2 = sci2(good2).feh
   mass2 = cat2(good2).logmass
   massmin2 = cat2(good2).logmass_min
   massmax2 = cat2(good2).logmass_max
   feherr2 = sci2(good2).feherr
   sn2 = sci2(good2).sn
   get_additional_uncertainties_lowmass,feh2,sn2,deltaage2,deltafeh_upper2,deltafeh_lower2,/sdss
   feherr_upper2 = sqrt(feherr2^2+deltafeh_upper2^2)
   feherr_lower2 = sqrt(feherr2^2+deltafeh_lower2^2)

   
   nit = 1000.
   mc_linfit,(mass-9),feh,feherr_upper,feherr_lower,nit,fitslope,fitintercept
   slope = fitslope[0]
   intercept = fitintercept[0]
   minslope = fitslope[0]-fitslope[1]
   maxslope = fitslope[0]+fitslope[1]
   minintercept = fitintercept[0]-fitintercept[1]
   maxintercept = fitintercept[0]+fitintercept[1]

   marr = findgen(16)/10.+8.
   enk_line = enk_slope*(marr-6)+enk_const
   nl_line = nl_slope*(marr-10)+nl_const
   best_line = slope[0]*(marr-9)+intercept[0]
   max_line = maxslope*(marr-9)+maxintercept
   min_line = minslope*(marr-9)+minintercept

   set_plot,'ps'
   !p.font=0
   sunsym = sunsymbol()
   psname='sdss_lowmass_mzr.eps'
   ymin = -1.2
   device,filename=psname,xsize=15,ysize=10,/encapsulated,/color
 ;     plot,mass,feh,psym=1,xtitle='Log(M/M'+sunsym+')',ytitle='[Fe/H]',/nodata,xrange=[8,9.5],yrange=[ymin,0.25],ystyle=1
      plot,mass,feh,psym=1,xtitle='log(mass)',ytitle='[Fe/H]',/nodata,xrange=[8,9.5],yrange=[ymin,0.25],ystyle=1,title='Measured with current SDSS data'
      polyfill,[marr,reverse(marr)],[(enk_line+enk_rms)>(ymin),(reverse(enk_line)-enk_rms)>(ymin)],color=fsc_color('blu2')
      oplot,marr,enk_line,linestyle=2,color=fsc_color('blu5')
      polyfill,[marr,reverse(marr)],[(nl_line+nl_rms)>(ymin),(reverse(nl_line)-nl_rms)>(ymin)],color=fsc_color('org2')
      oplot,marr,nl_line,linestyle=2,color=fsc_color('org5')

      cgerrplot,mass,feh-feherr_lower,feh+feherr_upper,color='darkgray'
      ;cgerrplot,feh,massmin,massmax,/horizontal
      ;oplot,mass(badmassmin),feh(badmassmin),psym=6
      ;oplot,mass(badmassmax),feh(badmassmax),psym=7
      vsym,4,rot=45,/fill
      oplot,mass,feh,psym=8,color=fsc_color('maroon')

      cgerrplot,mass2,feh2-feherr_lower2,feh2+feherr_upper2,color='darkgray'
      oplot,mass2,feh2,psym=cgsymcat(16),color=fsc_color('darkgreen'),symsize=0.8
      axis,yaxis=0,yrange=[ymin,0.25],ystyle=1,ytickformat='(A1)'
      axis,yaxis=1,yrange=[ymin,0.25],ystyle=1,ytickformat='(A1)'
      axis,xaxis=0,xrange=[8,9.5],xtickformat='(A1)' 
      al_legend,['15 < S/N per A < 25','10 < S/N per A < 15'],psym=[14,16],color=fsc_color(['maroon','darkgreen']),symsize=[1,0.8],box=0

;      oplot,marr,best_line,thick=2,color=fsc_Color('red')
;      oplot,marr,max_line,thick=2,color=fsc_Color('red'),linestyle=1
;      oplot,marr,min_line,thick=2,color=fsc_Color('red'),linestyle=1
   device,/close
   save,mass,massmin,massmax,feh,feherr_lower,feherr_upper,mass2,massmin2,massmax2,feh2,feherr_lower2,feherr_upper2,filename='lowmass_feh.sav'
   stop
end

