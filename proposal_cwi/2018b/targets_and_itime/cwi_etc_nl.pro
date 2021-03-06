function obj_cts,flux,texp
common assump, eff,skyflux,lambda,radius
   ;convert flux in erg to photon
  c_o = eff*(5.0341125d7*lambda*flux)*!dpi*radius^2*texp  ;photon/A
  return, c_o
end

function sky_cts,texp
common assump, eff,skyflux,lambda,radius
  c_s = eff*(5.0341125d7*lambda*skyflux)*!dpi*radius^2*texp  ;photon/A
  return, c_s
end

function cwi_etc_nl,flux,texp,snr_spatial_bin=snr_spatial_Bin
common assump, eff,skyflux,lambda,radius
;input: flux in erg/s/cm2/A/arcsec^2
;       snr_spatial_Bin in arcsec^2
;assume efficiency 10% source to detector
;assume sky brightness = 2.35e-17 erg/cm2/s/arcsec2/A
;assume lambda = 4500.
;assume radius = 5.1 meter
;assume skynoise dominated (ignore readnoise)
;return snr per angstrom

    eff = 0.1
    skyflux = 2.35d-17  ;erg/s/cm2/arcsec2/A
    lambda = 4500.
    radius = 5.1*100./2. ;cm
    read_noise = 2.5
    snr_spectral_bin = 1 
    if ~keyword_set(snr_spatial_bin) then snr_spatial_Bin = 1. 

    c_o = obj_cts(flux,texp)*snr_spatial_bin*snr_spectral_bin
    c_s = sky_cts(texp)*snr_spatial_bin*snr_spectral_bin
    ;pixel size is 15 micron per pixel, FOV is 60 arcsec per 24 mm detector
    ;1 pixel is 60*15e-3/24 = 0.0375 arcsec/pix
    pixels_per_arcsec2 = ((1./0.0375)^2)/4 ;2x2 binning
    c_r = read_noise*pixels_per_arcsec2*snr_spatial_bin


    snr = c_o/sqrt(c_s+c_o+c_r)
    ;stop
    return, snr

end
