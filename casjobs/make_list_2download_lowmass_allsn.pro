pro make_list_2download_lowmass_allsn
;This is to make proposal target list. (Need to measure HA EW of the galaxies)
cat=mrdfits('lowmass_allsn_nichaleet.fit',1)
good = where(cat.logmass_min gt 0. and cat.logmass_max gt 0. and cat.petromag_g lt 19,ngals)
print,'selecting ',ngals,'galaxies'
aa= cat(good)
plothist,aa.z
fname = 'list_lowmass_allsn_2downloadall.dat'
openw,1,fname
for i=0,ngals-1 do printf,1,[aa[i].plate,aa[i].mjd,aa[i].fiberid],format='(I6," ",I6," ",I6)'
close,1

stop
end
