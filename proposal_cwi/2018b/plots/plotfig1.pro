pro plotfig1
;GETTING SCI DATA
   sdss = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana_after_referee_report/sci_sdss_ana_afterev_final.fits',1)
   catsdss = mrdfits('/scr2/nichal/workspace2/ana/cl0024/newana/sci_sdss_ana.fits',2)

;restore values from leethochawalit 2018
   restore,'/scr2/nichal/workspace4/ana/ms0451/leetho18a_avevals.sav'
   sdss_slope = 0.16
   sdss_slopeerr = 0.02
   sdss_intercept = -0.05 ;at logm = 10
   sdss_intercepterr = 0.01
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;GETTING LITERATURE VALUES
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Get the average values measured in Gallazzi05
   g05_mass = [9.00,9.11,9.31,9.51,9.72,9.91,10.11,10.31,10.51,10.72,10.91,11.11,11.31,11.5];the first time is actually 8.91
   g05_feh  = [-0.6,-0.61,-0.65,-0.61,-0.52,-0.41,-0.23,-0.11,-0.01,0.04,0.07,0.10,0.12,0.13]
   ;below are the stars in the figure 8 of Gallazzi05
   g05_feherr = [0.62,0.56,0.59,0.55,0.47,0.43,0.35,0.31,0.27,0.25,0.22,0.21,0.2,0.2]/2.
   g05_fehlo = g05_feh-g05_feherr
   g05_fehhi = g05_feh+g05_feherr
   ;SDSS quiescent galaxies linear fit (Gallazzi14)
   g06_slope = 0.15
   g06_intercept = 0.109
   g06_rms = 0.14
   g06_mass = findgen(16)*0.1+10.
   g06_feh = (g06_mass-11.)*g06_slope+g06_intercept
   g06_Fehlo = g06_feh-g06_rms
   g06_fehhi = g06_feh+g06_rms

   toolow = where(g05_fehlo lt -1.,ctoolow)
   if ctoolow gt 0 then g05_fehlo(toolow) = -1

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Choi's Data
   z01 = {zlow:0.1,zhigh:0.2,mass:[9.9,10.2,10.4,10.7,11.0],Feh:[-0.05,-0.06,-0.01,-0.03,0.02],Feherr:[0.04,0.02,0.01,0.01,0.01]}
   z02 = {zlow:0.2,zhigh:0.3,mass:[10.2,10.5,10.7,11.0,11.3],Feh:[-0.08,-0.06,-0.03,-0.01,-0.05],Feherr:[0.04,0.02,0.01,0.01,0.02]}
   z03 = {zlow:0.3,zhigh:0.4,mass:[10.5,10.8,11.0,11.3],Feh:[-0.11,-0.05,-0.02,-0.03],Feherr:[0.03,0.01,0.01,0.02]}
   z04 = {zlow:0.4,zhigh:0.55,mass:[10.8,11.1,11.3],Feh:[-0.07,-0.04,-0.05],Feherr:[0.02,0.01,0.02]}
   z06 = {zlow:0.55,zhigh:0.7,mass:[10.9,11.0,11.3],Feh:[-0.15,-0.02,-0.05],Feherr:[0.07,0.03,0.03]}
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Sybilska et al 2017 (hELENa, IFU from Sauron data)
   aa=read_csv('/scr2/nichal/workspace2/catalogs/sybilska.csv',n_table_header=1,header=header)
   Syb_z0 = {mass:aa.field1,feh:aa.field2}
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;PLOTTING
   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   sunsym = sunsymbol()
   Delta = '!9'+string("104B)+'!x'
   alpha = '!9'+string("141B)+'!x'
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   psname='fig1_SDSS_FeH_mass.eps'
   device, filename = psname,xsize = 15,ysize = 10, $
   		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      ;make the outline of plots
      xrange=[9,11.5]
      yrange=[-1.,0.35]
      plot,sdss.logmstar,sdss.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
      
      ;shade the average region
      polyfill,[g06_mass,reverse(g06_mass)],[g06_fehhi,reverse(g06_fehlo)],color=fsc_color('pink')
      
      x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
      y=[hifeh_sdss,reverse(lofeh_sdss)]
      x(where(x eq 8.9)) = 9.
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=45
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=135
      rainbow_colors,n_colors=21
      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
      
      ;draw data points
      cgplot,sdss.logmstar,sdss.feh,psym=14,/overplot,color=40,symsize=1.3
      ;plot representative uncertainties
      lowermass = where(catsdss.logm50 lt 10.,complement=uppermass)
      cgerrplot,[-0.32,-0.32],[11.1+median(catsdss(lowermass).logm16-catsdss(lowermass).logm50),$
                               11.3+median(catsdss(uppermass).logm16-catsdss(uppermass).logm50)],$
                              [11.1+median(catsdss(lowermass).logm84-catsdss(lowermass).logm50),$
                               11.3+median(catsdss(uppermass).logm84-catsdss(uppermass).logm50)],$
                               /horizontal,color=40,thick=2
      cgerrplot,[11.1,11.3],[-0.32+median(sdss(lowermass).fehlower-sdss(lowermass).feh),$
                             -0.32+median(sdss(uppermass).fehlower-sdss(uppermass).feh)],$
                            [-0.32+median(sdss(lowermass).fehupper-sdss(lowermass).feh),$
                             -0.32+median(sdss(uppermass).fehupper-sdss(uppermass).feh)],$
                             color=40,thick=2

      
     ;Add Choi's data
      oploterror,z01.mass,z01.feh,z01.feherr,color=fsc_color('red5'),linethick=2,errcolor=10
      oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=fsc_color('red5'),symsize=2
      
      ;Add Sybilska2017
      oplot,syb_z0.mass,syb_z0.feh,color=fsc_color('maroon'),thick=3,linestyle=5
      
      ;Add best fitted line
      oplot,!x.crange,(!x.crange-10.)*sdss_slope+sdss_intercept,thick=3
       ;Labelling
      cglegend,title=['Leethochawalit et al. 2018','Gallazzi et al. 2014','Choi et al. 2014','Sybilska et al. 2017'],$
          psym=[14,15,46,0],location=[10.54,-0.65],box=0,charsize=0.8,/data,length=0,vspace=1.25,color=['navy','pink','red5','maroon'],symsize=1.2

      oplot,[10.34,10.48],[-0.65,-0.65],color=fsc_color('black'),thick=2
      oplot,[10.45,10.59],[-0.88,-0.88],linestyle=2,color=fsc_color('maroon'),thick=2
      oplot,[10.525],[-0.725],psym=cgsymcat(15),color=fsc_Color('pink'),symsize=1.3
      xyouts,9.1,0.2,'z~0 galaxies',charsize=1.
   device,/close

   psname='SDSS_FeH_mass_postdocapp.eps'
   !p.charsize=1.2
   device, filename = psname,xsize = 12,ysize = 8, $
   		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      ;make the outline of plots
      xrange=[9,11.5]
      yrange=[-1.,0.35]
      plot,sdss.logmstar,sdss.feh,/nodata,xrange=xrange,xstyle=5,yrange=yrange,ystyle=5
      
      ;shade the average region
      polyfill,[g06_mass,reverse(g06_mass)],[g06_fehhi,reverse(g06_fehlo)],color=fsc_color('pink')
      
      x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
      y=[hifeh_sdss,reverse(lofeh_sdss)]
      x(where(x eq 8.9)) = 9.
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=45
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=135
      rainbow_colors,n_colors=21
      ;draw axis
      axis,xaxis=0,xrange=xrange,xstyle=1,xtitle='Log(M/M'+sunsym+')'
      axis,yaxis=0,yrange=yrange,ystyle=1,ytitle='[Fe/H]'
      axis,xaxis=1,xrange=xrange,xstyle=1,xtickformat='(A1)'
      axis,yaxis=1,yrange=yrange,ystyle=1,ytickformat='(A1)'
      
      ;draw data points
      cgplot,sdss.logmstar,sdss.feh,psym=14,/overplot,color=40,symsize=1.3
      ;plot representative uncertainties
      lowermass = where(catsdss.logm50 lt 10.,complement=uppermass)
      cgerrplot,[-0.32,-0.32],[11.1+median(catsdss(lowermass).logm16-catsdss(lowermass).logm50),$
                               11.3+median(catsdss(uppermass).logm16-catsdss(uppermass).logm50)],$
                              [11.1+median(catsdss(lowermass).logm84-catsdss(lowermass).logm50),$
                               11.3+median(catsdss(uppermass).logm84-catsdss(uppermass).logm50)],$
                               /horizontal,color=40,thick=2
      cgerrplot,[11.1,11.3],[-0.32+median(sdss(lowermass).fehlower-sdss(lowermass).feh),$
                             -0.32+median(sdss(uppermass).fehlower-sdss(uppermass).feh)],$
                            [-0.32+median(sdss(lowermass).fehupper-sdss(lowermass).feh),$
                             -0.32+median(sdss(uppermass).fehupper-sdss(uppermass).feh)],$
                             color=40,thick=2

      
     ;Add Choi's data
      oploterror,z01.mass,z01.feh,z01.feherr,color=fsc_color('red5'),linethick=2,errcolor=10
      oplot,z01.mass,z01.feh,psym=cgsymcat(46),color=fsc_color('red5'),symsize=2
      
      ;Add Sybilska2017
      oplot,syb_z0.mass,syb_z0.feh,color=fsc_color('maroon'),thick=3,linestyle=5
      
      ;Add best fitted line
      oplot,!x.crange,(!x.crange-10.)*sdss_slope+sdss_intercept,thick=3
       ;Labelling
      cglegend,title=['Leethochawalit et al. 2018','Gallazzi et al. 2014','Choi et al. 2014','Sybilska et al. 2017'],$
          psym=[14,15,46,0],location=[10.04,-0.58],box=0,charsize=0.8,/data,length=0,vspace=1.25,$
          color=['navy','pink','red5','maroon'],symsize=[1.2]

      oplot,[9.84,9.98],[-0.58,-0.58],color=fsc_color('black'),thick=2
      oplot,[9.9,10.09],[-0.92,-0.92],linestyle=2,color=fsc_color('maroon'),thick=2
      xyouts,9.1,0.2,'z~0 galaxies',charsize=1.2
   device,/close

end
