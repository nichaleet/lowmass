function calsn_dbsp,mag,itime,lambda,grating,perang=perang
;calculate signal to noise per angstrom
;per pixel is the default
;mag is abmag
;lambda is in angstrom
;itime in second
    m_AB = mag   ;AB magnitude of object
    lambda = lambda * 1d-8    ;wavelength
    throughput = 0.12       ;total system throughput at wavelength of interest
    readnoise = 2.5 ;e per pixel
    if grating eq 600 then begin
        dlambda = 1.077d-8
        npix_dispersion = 1./1.077     ;number of pixels in 1 angstrom
    endif else if grating eq 1200 then begin
        dlambda = 4.10d-9
        npix_dispersion = 1./0.410     ;number of pixels in 1 angstrom
    endif

    texp = itime           ;exposure time in sec
    seeing_loss = 0.5      ;estimate of slit loss due to seeing
    npix_seeing = 2        ;about 2 pixels

    area = 1.76d5 ;cm^2
    h = 6.6260755d-27 ;erg s
    c = 2.99792458d10 ;cm s^-1
    
    f_lambda = 10.^(-0.4*(m_AB+48.6))/(h*lambda) ;nphotons/s/cm2/cm

    if keyword_set(perang) then begin
       photons = seeing_loss*texp*f_lambda*area*throughput;number of photons per cm
       photons = photons*1.d-8 ;number of photons per angstrom
       noise = sqrt(photons+readnoise*npix_seeing*npix_dispersion)
       snr = photons / noise 
    endif else begin 
       photons = texp*f_lambda*dlambda*area*throughput*seeing_loss ;nphotons/pixel
       noise = sqrt(photons+readnoise)
       snr = photons/noise
    endelse
return,snr 
end
