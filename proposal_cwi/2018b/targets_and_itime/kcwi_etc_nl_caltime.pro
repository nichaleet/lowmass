function obj_cts_kcwi,flux
common assump, eff,skyflux,lambda,radius
   ;convert flux in erg to photon
  c_o = eff*(5.0341125d7*lambda*flux)*!dpi*radius^2  ;photon/A/s
  return, c_o
end

function kcwi_etc_nl_caltime,flux,snr,snr_spatial_bin=snr_spatial_Bin,halfmoon=halfmoon
common assump, eff,skyflux,lambda,radius
;input: flux in erg/s/cm2/A/arcsec^2
;       snr_spatial_Bin in arcsec^2
;assume efficiency 10% source to detector
;assume sky brightness = 2.35e-17 erg/cm2/s/arcsec2/A
;assume lambda = 4500.
;assume radius = 5.1 meter
;assume skynoise dominated (ignore readnoise)
;return caltime

    eff = 0.2
    skyflux = 2.35d-17  ;erg/s/cm2/arcsec2/A https://lco.global/~apickles/Palomar/SkyBright/Aube.html
    if keyword_set(halfmoon) then skyflux = skyflux*6.3 ;2 mag brighter
    lambda = 4500.
    radius = 10*100./2. ;cm
    read_noise = 2.5
    snr_spectral_bin = 1 
    if ~keyword_set(snr_spatial_bin) then snr_spatial_Bin = 1. 

    s = obj_cts_kcwi(flux)*snr_spatial_bin*snr_spectral_bin

    n = eff*(5.0341125d7*lambda*skyflux)*!dpi*radius^2*snr_spatial_bin*snr_spectral_bin

    ;pixel size is 15 micron per pixel, FOV is 60 arcsec per 24 mm detector
    ;1 pixel is 60*15e-3/24 = 0.0375 arcsec/pix
    pixels_per_arcsec2 = ((1./0.0375)^2)/4 ;2x2 binning
    r = read_noise*pixels_per_arcsec2*snr_spatial_bin

    exptime = (snr^2*(s+n)+sqrt(SNR^4*(s+n)^2+4.*SNR^2*r*s^2))/(2*s^2)
    return, exptime

end
