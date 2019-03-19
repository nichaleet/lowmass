pro enk_nl_combine
  restore,'enk_table_template.sav'
  dwarf = read_ascii('enk_table.txt',template=enk_table_template) 
  restore,'/scr2/nichal/workspace2/ana/cl0024/newana/nl_sdss_feh.sav'
  restore,'/scr2/nichal/workspace2/lowmass/ana/lowmass_feh.sav'
 ;mass,massmin,massmax,feh,feherr_lower,feherr_upper,mass2,massmin2,massmax2,feh2,feherr_lower2,feherr_upper2

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
   marr = findgen(18)/2.+3
   enk_line = enk_slope*(marr-6)+enk_const
   nl_line = nl_slope*(marr-10)+nl_const


  set_plot,'ps'
  !p.font=0
  sunsym = sunsymbol()
  psname='enk_nl_mzr.eps'
  device,filename=psname,xsize=15,ysize=10,/encapsulated,/color
	plot,[dwarf.logm,NL_SDSS_MASS],[dwarf.feh,NL_SDSS_FEH],$
		/nodata,psym=1,xrange=[3,11.5],yrange=[-2.5,0.3],xstyle=1,ystyle=1
	polyfill,[8,9.5,9.5,8],[-2.5,-2.5,0.3,0.3],color=fsc_color('grn1')
        polyfill,[marr,reverse(marr)],[(enk_line+enk_rms)>(-2.5),(reverse(enk_line)-enk_rms)>(-2.5)],color=fsc_color('blu2')
        oplot,marr,enk_line,linestyle=2,color=fsc_color('blu5')
        polyfill,[marr,reverse(marr)],[(nl_line+nl_rms)>(-2.5),(reverse(nl_line)-nl_rms)>(-2.5)],color=fsc_color('org2')
        oplot,marr,nl_line,linestyle=2,color=fsc_color('org5')
	axis,xaxis=0,xtitle='Log(M/M'+sunsym+')',xrange=[3,11.5],xstyle=1
	axis,yaxis=0,ytitle='[Fe/H]',yrange=[-2.5,0.3],ystyle=1
        axis,yaxis=1,yrange=[-2.5,0.3],ystyle=1,ytickformat='(A1)'
	axis,xaxis=1,xrange=[3,11.5],xstyle=1,xtickformat='(A1)'
	oploterror,NL_SDSS_MASS,NL_SDSS_FEH,NL_SDSS_FEHERR,errcolor=fsc_color('gray'),psym=1
	cgerrplot,NL_SDSS_FEH,NL_SDSS_MASSLOW,NL_SDSS_MASSHIGH,/horizontal,color='gray'
	oploterror,dwarf.logm,dwarf.feh,dwarf.logm_err,dwarf.feh_err,errcolor=fsc_color('gray'),psym=1
	oplot,dwarf.logm,dwarf.feh,psym=cgsymcat(46),color=fsc_color('navy')
	oplot,NL_SDSS_MASS,NL_SDSS_FEH,psym=cgsymcat(46),color=fsc_color('org7')
   
        ;cgerrplot,mass,feh-feherr_lower,feh+feherr_upper,color='gray'
        ;oplot,mass,feh,psym=cgsymcat(14),color=fsc_color('maroon'),symsize=0.5
        ;cgerrplot,mass2,feh2-feherr_lower2,feh2+feherr_upper2,color='gray'
        ;oplot,mass2,feh2,psym=cgsymcat(16),color=fsc_color('darkgreen'),symsize=0.5

	al_legend,['Kirby et al. 2013','Leethochawalit et al. (in prep)'],psym=46,color=fsc_color(['navy','org7']),/right,/bottom,box=0
  device,/close
  stop	 
end
