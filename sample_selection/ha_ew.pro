pro gaussfunc,x,a,f,pder
   z = (x-a[1])/a[2]
   f = a[0]*exp(-z^2/2.)+a[3]
   pder = [[exp(-z^2/2.)],[a[0]*exp(-z^2/2.)*(x-a[1])/(a[2])^2],[a[0]*exp(-z^2/2.)*z^2/a[2]],[replicate(1.,n_elements(x))]]
end

function ha_ew,flux,lambda,ivar,z,ewerr,nostop=nostop
   fixz=0
   lambdaobs = lambda
   while fixz eq 0 do begin
      lamrange=[6553.,6580.]
      lambdarest  =lambdaobs/(z+1.)
      inrange = where(lambdarest gt lamrange[0] and lambdarest le lamrange[1] and finite(flux) and finite(ivar),cinrange)
      floorlev = median(flux(where(lambdarest gt 6520. and lambdarest lt 6620)))
      dummy = min(abs(lambdarest-6564.6),locpeak)
      peak = flux(locpeak)-floorlev
      coeff = [peak,6564.6,1.,floorlev] 
      yfit = curvefit(lambdarest(inrange),flux(inrange),ivar(inrange),coeff,sigmacoeff,fita=[1,0,1,1],$
                      function_name='gaussfunc',status=status)
      if status eq 0 then print,'fit success' else print, 'fit failed'
      print, coeff
      EW = sqrt(2.*!pi)*coeff[0]*abs(coeff[2])/coeff[3]
      EWerr = abs(EW)*sqrt((sigmacoeff[0]/coeff[0])^2+(sigmacoeff[2]/coeff[2])^2+(sigmacoeff[3]/coeff[3])^2)
      	
      ;plotting
      rangeplot=[6530,6610]
      plot,lambdarest,flux,xrange=rangeplot,/nodata
      oplot,lambdarest,flux,psym=10
      oplot,lambdarest(inrange),yfit,color=fsc_color('yellow')
      oplot,[6564,6564],!y.crange,color=fsc_Color('red'),linestyle=1
      xyouts, 0.2,0.3,'EW= '+sigfig(EW,3)+'+/-'+sigfig(EWerr,3),/normal
      if ~keyword_set(nostop) then begin
         strfixz=''
         read,strfixz,prompt='fix redshift? (y/n)'
         if strfixz eq 'y' then fixz=0 else fixz=1
         if fixz eq 0 then begin
            read,newhawl,prompt='type in the location of Ha that you see: '
      	    z = (newhawl*(z+1)/6564.)-1.
         endif
       endif else fixz=1
   endwhile
   return,EW
end
