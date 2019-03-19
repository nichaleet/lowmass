pro enk_nl_combine
  restore,'/scr2/nichal/workspace_lowmass/proposal_dbsp/enk_table_template.sav'
  dwarf = read_ascii('/scr2/nichal/workspace_lowmass/proposal_dbsp/enk_table.txt',template=enk_table_template) 

  sdss = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/sci_sdss_ana_afterev_final.fits',1)
  catsdss = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana/sci_sdss_ana.fits',2)
                                                       
   enk_slope = 0.3   ;dex                              
   enk_slope_err = 0.02 ;dex                           
   enk_rms = 0.17    ;dex                              
   enk_const = -1.69 ;dex this is constant at 10^6 Msu n
   enk_const_err = 0.04 ;dex
   nl_rms = 0.13
   nl_slope = 0.16   ;dex
   nl_slope_err = 0.03 ;dex
   nl_const = -0.05 ;dex this is constant at 10^10 msun (y=m(x-10)+c)
   nl_const_err = 0.01 ;dex
   marr = findgen(18)/2.+3
   enk_line = enk_slope*(marr-6)+enk_const
   nl_line = nl_slope*(marr-10)+nl_const


  set_plot,'ps'
  !p.font=0
  sunsym = sunsymbol()
  psname='enk_nl_mzr.eps'
  device,filename=psname,xsize=15,ysize=10,/encapsulated,/color
	plot,[dwarf.logm,sdss.logmstar],[dwarf.feh,sdss.feh],$
		/nodata,psym=1,xrange=[3,11.5],yrange=[-2.5,0.3],xstyle=1,ystyle=1
	polyfill,[8,9.,9.,8],[-2.5,-2.5,0.3,0.3],color=fsc_color('grn1')
        polyfill,[marr,reverse(marr)],[(enk_line+enk_rms)>(-2.5),(reverse(enk_line)-enk_rms)>(-2.5)],color=fsc_color('org2')
        oplot,marr,enk_line,linestyle=2,color=fsc_color('org5')
        polyfill,[marr,reverse(marr)],[(nl_line+nl_rms)>(-2.5),(reverse(nl_line)-nl_rms)>(-2.5)],color=fsc_color('blu2')
        oplot,marr,nl_line,linestyle=2,color=fsc_color('blu5')
	axis,xaxis=0,xtitle='Log(M/M'+sunsym+')',xrange=[3,11.5],xstyle=1
	axis,yaxis=0,ytitle='[Fe/H]',yrange=[-2.5,0.3],ystyle=1
        axis,yaxis=1,yrange=[-2.5,0.3],ystyle=1,ytickformat='(A1)'
	axis,xaxis=1,xrange=[3,11.5],xstyle=1,xtickformat='(A1)'
        cgerrplot,sdss.logmstar,sdss.fehlower,sdss.fehupper,color='gray'
	cgerrplot,sdss.feh,catsdss.logm16,catsdss.logm84,/horizontal,color='gray'
	oploterror,dwarf.logm,dwarf.feh,dwarf.logm_err,dwarf.feh_err,errcolor=fsc_color('gray'),psym=1
	oplot,dwarf.logm,dwarf.feh,psym=cgsymcat(46),color=fsc_color('org7')
	oplot,sdss.logmstar,sdss.feh,psym=cgsymcat(46),color=fsc_color('navy')
   
        ;cgerrplot,mass,feh-feherr_lower,feh+feherr_upper,color='gray'
        ;oplot,mass,feh,psym=cgsymcat(14),color=fsc_color('maroon'),symsize=0.5
        ;cgerrplot,mass2,feh2-feherr_lower2,feh2+feherr_upper2,color='gray'
        ;oplot,mass2,feh2,psym=cgsymcat(16),color=fsc_color('darkgreen'),symsize=0.5

	al_legend,['Kirby et al. 2013','Leethochawalit et al. 2018'],psym=46,color=fsc_color(['navy','org7']),/right,/bottom,box=0
  device,/close

  ;the following is for postdoc app
   ;lowmass stuff
   sci = mrdfits('/scr2/nichal/workspace_lowmass/sps_fit/data/sdss_lowmass/sps_fit03.fits.gz',1) ;this spec was smoothed with v=200 km/s so S/N is off
   sciforsn = mrdfits('/scr2/nichal/workspace_lowmass/sps_fit/data/sdss_lowmass/sps_fit01.fits.gz',1)
   cat = mrdfits('/scr2/nichal/workspace_lowmass/proposal_cwi/targets_and_itime/combine_targets.fits',1)
   cat2 = mrdfits('/scr2/nichal/workspace_lowmass/proposal_cwi/targets_and_itime/combine_targets.fits',2)
   symsize = 0.5*(sciforsn.sn/(max(sciforsn.sn)-min(sciforsn.sn)))+0.5
   radfraction = (3./cat2.petrorad_g)  ;sdss fiber is 3 arcsec
   color=bytscl(radfraction,min=0,max=0.5)
   feherr_upper = sqrt(sci.feherr^2+(sci.fehupper-sci.feh)^2)
   feherr_lower = sqrt(sci.feherr^2+(sci.fehlower-sci.feh)^2)


  !p.charsize=1.2
  psname='enk_nl_mzr_postdocapp.eps'
  device,filename=psname,xsize=12,ysize=8,/encapsulated,/color
	plot,[dwarf.logm,sdss.logmstar],[dwarf.feh,sdss.feh],$
		/nodata,psym=1,xrange=[3,11.5],yrange=[-2.5,0.3],xstyle=1,ystyle=1
	polyfill,[8,9.,9.,8],[-2.5,-2.5,0.3,0.3],color=fsc_color('grn1')
        polyfill,[marr,reverse(marr)],[(enk_line+enk_rms)>(-2.5),(reverse(enk_line)-enk_rms)>(-2.5)],color=fsc_color('org2')
        oplot,marr,enk_line,linestyle=2,color=fsc_color('org5')
        polyfill,[marr,reverse(marr)],[(nl_line+nl_rms)>(-2.5),(reverse(nl_line)-nl_rms)>(-2.5)],color=fsc_color('blu2')
        oplot,marr,nl_line,linestyle=2,color=fsc_color('blu5')
	axis,xaxis=0,xtitle='Log(M/M'+sunsym+')',xrange=[3,11.5],xstyle=1
	axis,yaxis=0,ytitle='[Fe/H]',yrange=[-2.5,0.3],ystyle=1
        axis,yaxis=1,yrange=[-2.5,0.3],ystyle=1,ytickformat='(A1)'
	axis,xaxis=1,xrange=[3,11.5],xstyle=1,xtickformat='(A1)'
        cgerrplot,sdss.logmstar,sdss.fehlower,sdss.fehupper,color='gray'
	cgerrplot,sdss.feh,catsdss.logm16,catsdss.logm84,/horizontal,color='gray'
	oploterror,dwarf.logm,dwarf.feh,dwarf.logm_err,dwarf.feh_err,errcolor=fsc_color('gray'),psym=1
	oplot,dwarf.logm,dwarf.feh,psym=cgsymcat(46),color=fsc_color('org7')
	oplot,sdss.logmstar,sdss.feh,psym=cgsymcat(46),color=fsc_color('navy')
        cgerrplot,sci.logmstar,sci.feh-feherr_lower,sci.feh+feherr_upper,color='darkgray'
        for i=0,n_elements(sci)-1 do oplot,[sci[i].logmstar],[sci[i].feh],psym=cgsymcat(46),symsize=symsize[i],color=fsc_color('maroon')
   
        ;cgerrplot,mass,feh-feherr_lower,feh+feherr_upper,color='gray'
        ;oplot,mass,feh,psym=cgsymcat(14),color=fsc_color('maroon'),symsize=0.5
        ;cgerrplot,mass2,feh2-feherr_lower2,feh2+feherr_upper2,color='gray'
        ;oplot,mass2,feh2,psym=cgsymcat(16),color=fsc_color('darkgreen'),symsize=0.5
        xyouts,3.5,0.,'z~0 galaxies'
	al_legend,['Kirby et al. 2013','Leethochawalit et al. 2018','From current SDSS spectra'],psym=46,color=['org7','navy','maroon'],position=[6.1,-1.8],box=0,charsize=0.9
  device,/close


  stop	 
end

