pro sagaana
   cat = mrdfits('saga_spectra_clean.fits',1,hdr)  
   good = where(aa.spec_ha_ew gt -90 and aa.spec_ha_ew lt 0.,cgood)
   ;calculate mass 
end
