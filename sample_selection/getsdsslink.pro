pro getsdsslink
  str = mrdfits('lowmass_targets.fits',1)
  for i=0,n_elements(str)-1 do begin
     aa=mrdfits(str[i].specfile,2,/silent)
     print,'http://skyserver.sdss.org/dr12/en/tools/quicklook/summary.aspx?id='+aa.bestobjid
  endfor
end
