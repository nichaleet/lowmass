pro make_worm_plot_lowmass,dir

	fehplot = [-0.8,-0.6,-0.4,-0.2]
	ageplot = [2.5,3.5,4.5,5.5,6.5,7.5]
	nfeh = n_Elements(fehplot)
	
	set_plot,'ps'
	psname=dir+'/chisq_worm_plot_'+dir+'.eps'
	!p.font=0
	device, filename = psname,xsize = 15,ysize = 10, $
		xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
	plot,findgen(601)/500.-1.2,alog10(findgen(12)+0.5)+9.,xrange=[-1.2,0.],yrange=alog10([2,12.])+9,/nodata,$
		xtitle='[Fe/H]',xstyle=1,ystyle=5;,position=[0.15,0.1,0.85,0.9]
	axis, yaxis=0,ytitle='log(Age/yr)',yrange=alog10([2,12.])+9,ystyle=1
	axis, yaxis=1,yticks=6,ytickv=alog10([2,3,4,5,6,8,10.])+9.,ytickn=['2','3','4','5','6','8','10'],ycharsize=0.7
	burd_colors
	fehcolor=[10,60,190,254]
	for i=0,nfeh-1 do for j=0,n_elements(ageplot)-1 do begin
		fehref = fehplot(i)
		ageref = ageplot(j)
		;read data
		fitsfile = dir+"/chisq_age"+sigfig(ageref,2)+"_feh"+sigfig(fehref,4)+".fits"	
		if file_test(fitsfile) eq 0 then stop,'no fits file found'
		chisqarr = mrdfits(fitsfile,0,header)
		;you should check below if it matches those in degeneracy.pro
		grid_feh = rebin(findgen(601)/500.-1.2,601,12) ;range from [-0.6,0.192] dex
		grid_age = transpose(rebin(alog10(findgen(12)+0.5)+9.,12,601)) ;range from [0.5,13] Gyr
		minchisq = min(chisqarr)
		contour,chisqarr,grid_feh,grid_age,levels=[minchisq+0.02,minchisq+1],c_linestyle=0,/overplot,color=fehcolor(i)
		oplot,[fehref],alog10([ageref])+9.,psym=7,color=fehcolor(i)
		
	endfor
	device,/close
	stop
end
