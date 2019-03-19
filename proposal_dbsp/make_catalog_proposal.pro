pro make_catalog_proposal,ra,dec,logmstar,gmag,sn
    openw, lun2, 'target_list.tex', /get_lun
    printf, lun2, '\begin{deluxetable*}{llllll}'
    printf, lun2, '\tablewidth{0pt}'
    printf, lun2, '\tablecolumns{6}'
    printf, lun2, '\tablehead{\colhead{No.} & \colhead{RA} & \colhead{DEC} & \colhead{$log(M_*/M_\odot)$} & \colhead{g Magnitude} & \colhead{current SDSS SN}'

    printf, lun2, '\startdata'
    
    n = n_elements(ra)
    ;only works when n=72
    for i=0, n-1 do begin
        mags =  sigfig(gmag[i],3)
        radec, ra[i], dec[i], r1, r2, r3, d1, d2, d3
        ras = string(r1, r2, r3, format='(I02,1X,I02,1X,D05.2)')
        decs = (dec[i] lt 0 ? '-' : '+')+string(abs(d1), abs(d2), abs(d3), format='(I02,1X,I02,1X,D04.1)')
	mass = string(logmstar[i],format='(D4.1)')
        if i ge 99 then nos = string(i+1,format='(I03)')
        if i ge 9 and i lt 99 then nos = string(i+1,format='(I02)') 
        if i lt 9 then nos= string(i+1,format='(I01)')
       ; printf, lun2, nos, ras, decs, mass, mags,sn[i], format='(A3,2(" & ",A11)," & ",A5," & ",A5," & ",D5.1," \\")'
 printf, lun2, nos, ras, decs, mass, mags,sn[i], format='(A3,2(" & ",A11)," & ",A5," & ",A5," \\")'
    endfor
    close,lun2
    free_lun, lun2
end
