function calitime_dbsp,mag,sn,lambda,grating
;calculate signal to noise per angstrom
;SN to noise per angstrom is the default
;mag is abmag
;lambda is in angstrom
;itime in second
    m_AB = mag   ;AB magnitude of object
    lambda = lambda * 1d-8    ;wavelength
    throughput = 0.12       ;total system throughput at wavelength of interest
    readnoise = 2.5 ;e per pixel
    skymag = 20.9 ;mag per arcsec^2 in gband https://lco.global/~apickles/Palomar/SkyBright/Aube.html

    if grating eq 600 then begin
        dlambda = 1.077d-8
        npix_dispersion = 1./1.077     ;number of pixels in 1 angstrom
    endif else if grating eq 1200 then begin
        dlambda = 4.10d-9
        npix_dispersion = 1./0.410     ;number of pixels in 1 angstrom
    endif

    seeing_loss = 0.5      ;estimate of slit loss due to seeing
    npix_seeing = 2        ;about 2 pixels

    area = 1.76d5 ;cm^2
    h = 6.6260755d-27 ;erg s
    c = 2.99792458d10 ;cm s^-1
    
    f_lambda = 10.^(-0.4*(m_AB+48.6))/(h*lambda) ;nphotons/s/cm2/cm
    f_sky = 10.^(-0.4*(skymag+48.6))/(h*lambda) ;nphotons/s/cm2/cm/arcsec^2

   pps = seeing_loss*f_lambda*area*throughput;number of photons per cm per second
   pps = pps*1.d-8 ;number of photons per angstrom per second
   sps = f_sky*area*throughput*1.d-8 ;sky photons per second per angstrom assum 1 arcsec^2
   
   RN  = readnoise*npix_seeing*npix_dispersion ;readnoise per read
   itime = ((SN^2*(pps+sps))+sqrt((SN^4*(pps+sps)^2)+(4.*pps^2*RN*SN^2)))/(2.*pps^2)
   
return,itime
end
