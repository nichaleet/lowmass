pro make_targetlist
   list = mrdfits('lowmass_targets.fits',1)
   radec,list.ra,list.dec,ihr,imin,xsec,ideg,imn,xsc
   nobjs = n_elements(ihr)


   openw, lun2, 'leethochawalit_dbsp_180513_targetlist.txt', /get_lun
   pa = [102,0,65,26,176,89,124,170,25,103,171,21,11,162,136,72,175,29,23,13,0,61,81,20,135,132,130,349,120,72,68,$
         101,159,98,175,187,5,8,45,28,68,118,102,91,79,57,158,150,68,0,45,162,32,61,120,94,127,176,68,148,120,15,7,$
         82,27,140,6,96,90,84,0,125,126]
   for i=0,nobjs-1 do begin
      ra = string(ihr[i],imin[i],xsec[i], format='(I02,1X,I02,1X,D05.2)')
      dec = (list[i].dec lt 0 ? '-' : '+')+string(abs(ideg[i]), abs(imn[i]), abs(xsc[i]), format='(I02,1X,I02,1X,D04.1)')
      name =  string(i+1,format='("lowmass",I03)')  
      printf,lun2,name,ra,dec,pa[i],format='(A,",",A,",",A,",2000,",I0)'
   endfor

   list = mrdfits('highermass_targets.fits',1)
   radec,list.ra,list.dec,ihr,imin,xsec,ideg,imn,xsc
   nobjs = n_Elements(ihr)
   for i=0,nobjs-1 do begin
      ra = string(ihr[i],imin[i],xsec[i], format='(I02,1X,I02,1X,D05.2)')
      dec = (list[i].dec lt 0 ? '-' : '+')+string(abs(ideg[i]), abs(imn[i]), abs(xsc[i]), format='(I02,1X,I02,1X,D04.1)')
      name =  string(i+74,format='("lowmass",I03)')  
      printf,lun2,name,ra,dec,format='(A,",",A,",",A,",2000")'
   endfor

   close,lun2
   free_lun, lun2

end
