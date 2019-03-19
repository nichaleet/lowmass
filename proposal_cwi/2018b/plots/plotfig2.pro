pro plotfig2

   enk_slope = 0.3   ;dex
   enk_slope_err = 0.02 ;dex
   enk_rms = 0.17    ;dex
   enk_const = -1.69 ;dex this is constant at 10^6 Msun
   enk_const_err = 0.04 ;dex
   nl_rms = 0.13
   nl_slope = 0.16   ;dex
   nl_slope_err = 0.02 ;dex
   nl_const = -0.05 ;dex this is constant at 10^10 msun (y=m(x-10)+c)
   nl_const_err = 0.01 ;dex

   ;lowmass stuff
   sci = mrdfits('/scr2/nichal/workspace_lowmass/sps_fit/data/sdss_lowmass/sps_fit03.fits.gz',1) ;this spec was smoothed with v=200 km/s so S/N is off
   sciforsn = mrdfits('/scr2/nichal/workspace_lowmass/sps_fit/data/sdss_lowmass/sps_fit01.fits.gz',1)
   cat = mrdfits('/scr2/nichal/workspace_lowmass/proposal_cwi/targets_and_itime/combine_targets.fits',1)
   cat2 = mrdfits('/scr2/nichal/workspace_lowmass/proposal_cwi/targets_and_itime/combine_targets.fits',2)
   symsize = 1.5*(sciforsn.sn/(max(sciforsn.sn)-min(sciforsn.sn)))+1.
   radfraction = (3./cat2.petrorad_g)  ;sdss fiber is 3 arcsec
   color=bytscl(radfraction,min=0,max=0.5)
   feherr_upper = sqrt(sci.feherr^2+(sci.fehupper-sci.feh)^2)
   feherr_lower = sqrt(sci.feherr^2+(sci.fehlower-sci.feh)^2)
  
;   nit = 1000.
;   mc_linfit,(mass-9),feh,feherr_upper,feherr_lower,nit,fitslope,fitintercept
;   slope = fitslope[0]
;   intercept = fitintercept[0]
;   minslope = fitslope[0]-fitslope[1]
;   maxslope = fitslope[0]+fitslope[1]
;   minintercept = fitintercept[0]-fitintercept[1]
;   maxintercept = fitintercept[0]+fitintercept[1]

   marr = findgen(11)/10.+8.
   enk_line = enk_slope*(marr-6)+enk_const
   nl_line = nl_slope*(marr-10)+nl_const
;   best_line = slope[0]*(marr-9)+intercept[0]
;   max_line = maxslope*(marr-9)+maxintercept
;   min_line = minslope*(marr-9)+minintercept

   set_plot,'ps'
   !p.font=0
   sunsym = sunsymbol()
   psname='fig2_sdss_lowmass_mzr.eps'
   ymin = -1.2
   device,filename=psname,xsize=15,ysize=10,/encapsulated,/color
      plot,sci.logmstar,sci.feh,psym=1,xtitle='Log(M/M'+sunsym+')',ytitle='[Fe/H]',/nodata,xrange=[8,9.],yrange=[ymin,0.25],ystyle=1,title='Measured with current SDSS data'
      polyfill,[marr,reverse(marr)],[(enk_line+enk_rms)>(ymin),(reverse(enk_line)-enk_rms)>(ymin)],color=fsc_color('org2')
      oplot,marr,enk_line,linestyle=2,color=fsc_color('org5')
      polyfill,[marr,reverse(marr)],[(nl_line+nl_rms)>(ymin),(reverse(nl_line)-nl_rms)>(ymin)],color=fsc_color('blu2')
      oplot,marr,nl_line,linestyle=2,color=fsc_color('blu5')

      cgerrplot,sci.logmstar,sci.feh-feherr_lower,sci.feh+feherr_upper,color='darkgray'
;      cgerrplot,sci.feh,sci.logm16,sci.logm84,/horizontal,color='darkgray'
      ;oplot,mass(badmassmin),feh(badmassmin),psym=6
      ;oplot,mass(badmassmax),feh(badmassmax),psym=7
      ;cgloadct, 33
      ;for i=0,n_elements(sci)-1 do oplot,[sci[i].logmstar],[sci[i].feh],psym=cgsymcat(46),symsize=symsize[i],color=color[i]
      for i=0,n_elements(sci)-1 do oplot,[sci[i].logmstar],[sci[i].feh],psym=cgsymcat(46),symsize=symsize[i],color=fsc_color('maroon')

      axis,yaxis=0,yrange=[ymin,0.25],ystyle=1,ytickformat='(A1)'
      axis,yaxis=1,yrange=[ymin,0.25],ystyle=1,ytickformat='(A1)'
      axis,xaxis=0,xrange=[8,9.],xtickformat='(A1)' 
;      al_legend,['15 < S/N per A < 25','10 < S/N per A < 15'],psym=[14,16],color=fsc_color(['maroon','darkgreen']),symsize=[1,0.8],box=0

;      oplot,marr,best_line,thick=2,color=fsc_Color('red')
;      oplot,marr,max_line,thick=2,color=fsc_Color('red'),linestyle=1
;      oplot,marr,min_line,thick=2,color=fsc_Color('red'),linestyle=1
   device,/close

   device,filename='fig2_rpetrohist.epx',xsize=15,ysize=10,/encapsulated,/color
      !p.font=0
      plothist,cat2.petrorad_g,xtitle='Petrosian Radius (arcsec)',xrange=[0,30],bin=2
   device,/close
   ;calculate effective radii base on Graham2005

   stop
end

