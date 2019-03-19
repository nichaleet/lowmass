pro mc_linfit,x,y,yerr_up,yerr_low,nit,slope,intercept
 ;monte carlo linear fit 
 slopearr = fltarr(nit)
 constarr = fltarr(nit)
 nnum = n_elements(x)
 meanerr = 0.5*(yerr_low+yerr_up)
 for i=0,nit-1 do begin
    ysamp = randomu(seed,nnum)*(yerr_up+yerr_low)+(y-yerr_low)
    param = linfit(x,ysamp);,measure_Errors=meanerr)
    slopearr(i) = param[1]
    constarr(i) = param[0]
 endfor
 slope = [mean(slopearr),stdev(slopearr)]
 intercept = [mean(constarr),stdev(constarr)]
end
