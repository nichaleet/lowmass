pro make_list_2download
aa=mrdfits('lowmass_nichaleet.fit',1)
pothist,aa.z
ngals = n_elements(aa)
fname = 'list_lowmass_2download.dat'
openw,1,fname
for i=0,ngals-1 do printf,1,[aa[i].plate,aa[i].mjd,aa[i].fiberid],format='(I6," ",I6," ",I6)'
close,1


end
