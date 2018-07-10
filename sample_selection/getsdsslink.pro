pro getsdsslink
  str = mrdfits('highermass_targets.fits',1)
  for i=0,n_elements(str)-1 do begin
     aa=mrdfits(str[i].specfile,2,/silent)
     if tag_exist(aa,'bestobjid') then $
         print,'http://skyserver.sdss.org/dr12/en/tools/quicklook/summary.aspx?id='+aa.bestobjid else $
         print,'http://skyserver.sdss.org/dr12/en/tools/quicklook/summary.aspx?id='+aa.objid
  endfor

end
