pro ngc147itime

   ;use ds9 to measure counts in 10''x40'' box in hst_06233_02_wfpc2_f555w_wf_drz.fits
   ;calculate magniftude and surface brightness from the counts and zeropoint info in the header
   ;See the regions in obj_cwiplan.reg (planned for 3 pointings: north, obj, south)

   clight = 2.99792458d10 ;cgs
   photzpt = -21.1   ;from header
   photflam = 4.891731741573034E-19 ;from header
   zeropt = -2.5*alog10(photflam)+photzpt 
   lambda = 5443.
   snrarr = [10.,20.,30.,50.]
   area = 400. ;arcsec2
   obj = replicate({name:'',mag:0.,magerr:0.,area:area,sfb:0.,sfberr:0.,sff:0.,sfferr:0.,$
                        snr:snrarr,itime:fltarr(n_Elements(snrarr))},9)
   obj.name = ['c1','c2','c3','c4','c5','c6','n6','s1','s3']
   obj.area = [area,area,area,area,area,area,area*2.,area*2.,area*2.]
   ccount = [2271,2709,3442,6882,2964,2443,1221,1403,610]
   ccounterr = [48,52,59,83,54,49,35,38,25]
   csub = [0,0,876,3510+500,500,0.,0.,0.,0.]
   csuberr = [0,0,30,sqrt(60^2+15^2),15,0.,0.,0.,0.]
   ccount = ccount-csub
   ccounterr = sqrt(ccounterr^2+csuberr^2)


   ;calculate magnitude of each slice
   obj.mag = -2.5*alog10(ccount)+zeropt
   obj.magerr = abs(2.5*(ccounterr/ccount/alog(10.)))

   flux = 10.^((obj.mag+21.1)/(-2.5)) ;erg/cm2/s/ang
   fluxerr = abs(flux*alog(10.)*obj.magerr/2.5)
   
   obj.sff = flux/area ;erg/cm2/s/ang/arcsec2
   obj.sfferr = fluxerr/area

   obj.sfb = -2.5*alog10(obj.sff)-21.1
   obj.sfberr = abs(2.5*(obj.sfferr/obj.sff/alog(10.)))

   for i=0,n_elements(snrarr)-1 do begin
       snr = snrarr[i]
       for j=0,n_elements(obj)-1 do begin
          sff = obj[j].sff
          itime = cwi_etc_nl_caltime(sff,snr,snr_spatial_bin=obj[j].area)    
          obj[j].itime[i] =itime
       endfor
   endfor
   print, 'snr   ', snrarr
   for i=0,n_elements(obj)-1 do begin
      print, obj[i].name,obj[i].itime
   endfor
stop
end
