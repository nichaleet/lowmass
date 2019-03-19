pro mkcwiregion,ra,dec,reff,r50,name
   rastring = strtrim(string(ra,format='(f11.7)'),2)
   decstring = strtrim(string(dec,format='(f11.8)'),2)
   openw,1,name
   printf,1,'# Region file format: DS9 version 4.1'
   printf,1,'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
   printf,1,'fk5'
   printf,1,'circle('+rastring+','+decstring+','+strtrim(string(reff,format='(f6.3)'),2)+'")'
   printf,1,'circle('+rastring+','+decstring+','+strtrim(string(r50,format='(f6.3)'),2)+'") # color=magenta dash=1'
   printf,1,'box('+rastring+','+decstring+',60.000",40.000",0) # color=blue width=2'
   printf,1,'circle('+rastring+','+decstring+',3.000") # color=red'
   close,1
end

pro fig4_createregion
   cat = mrdfits('/scr2/nichal/workspace_lowmass/proposal_cwi/targets_and_itime/combine_targets.fits',1)
   cat2 = mrdfits('/scr2/nichal/workspace_lowmass/proposal_cwi/targets_and_itime/combine_targets.fits',2)
   cat3 = mrdfits('/scr2/nichal/workspace_lowmass/proposal_cwi/targets_and_itime/combine_targets.fits',3)
   ;create region for all objects
   directory = '/scr2/nichal/workspace_lowmass/img/'
   for i=0,n_elements(cat)-1 do begin
       name = directory+'obj'+strtrim(string(i),2)+'-'+strtrim(string(cat[i].plate),2)+'-'+$
              strtrim(string(cat[i].mjd),2)+'-'+strtrim(string(cat[i].fiberid),2)+'.reg'
       mkcwiregion,cat[i].ra,cat[i].dec,cat2[i].petror90_g,cat2[i].petror50_g,name
   endfor


end
