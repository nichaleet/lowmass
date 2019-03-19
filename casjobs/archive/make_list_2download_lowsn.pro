pro make_list_2download_lowsn
cat=mrdfits('lowmass_lowsn_nichaleet.fit',1)
good = where(cat.logmass_min gt 0. and cat.logmass_max gt 0. and cat.petromag_g lt 19 and cat.logmass lt 8.8 and cat.snmedian gt 10.,ngals)
print,'selecting ',ngals,'galaxies'
aa= cat(good)
plothist,aa.z
fname = 'list_lowmass_lowsn_2download.dat'
openw,1,fname
for i=0,ngals-1 do printf,1,[aa[i].plate,aa[i].mjd,aa[i].fiberid],format='(I6," ",I6," ",I6)'
close,1


end
