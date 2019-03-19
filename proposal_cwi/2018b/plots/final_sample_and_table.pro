pro final_sample_and_table
   cat = mrdfits('/scr2/nichal/workspace_lowmass/proposal_cwi/targets_and_itime/combine_targets.fits',1)
   cat2 = mrdfits('/scr2/nichal/workspace_lowmass/proposal_cwi/targets_and_itime/combine_targets.fits',2)
   cat3 = mrdfits('/scr2/nichal/workspace_lowmass/proposal_cwi/targets_and_itime/combine_targets.fits',3)

   good = where(cat3.PETRO90ITIME lt 3600. and cat2.PETROR90_G gt 5. and cat2.PETROR90_G lt 25.,cgood)
   good = good(where(good lt 82 or good gt 84 and good ne 41))
   cgood = n_Elements(good)
   cat = cat(good)
   cat2 = cat2(good)
   cat3 = cat3(good)

   factor = fltarr(cgood)+2.
   wlarge = where(cat2.petror90_g gt 20.)
   factor(wlarge) = 4.
   total_itime = cat3.petro90itime*factor
   all_itime = total(total_itime)/3600.
   openw, lun2, 'cwi_target_list.tex', /get_lun
   printf, lun2, '\begin{tabular}{lllllllll}'
   printf, lun2, 'No. & RA & DEC & g Mag& $log(M_*/M_\odot)$ & flux$_{50}$ & flux$_{90}$ & R$_{\text{eff}}$ & Itime(s) & Total Itime(s) \\'
   printf, lun2, '\hline'

   for i=0,cgood-1 do begin
      radec, cat[i].ra, cat[i].dec, r1, r2, r3, d1, d2, d3
      ras = string(r1, r2, r3, format='(I02,1X,I02,1X,D05.2)')
      decs = (cat[i].dec lt 0 ? '-' : '+')+string(abs(d1), abs(d2), abs(d3), format='(I02,1X,I02,1X,D04.1)')
      mass = string(cat[i].logmass,format='(D4.2)')
      gmag = string(cat[i].petromag_g,format='(D4.1)')
      flux50 = string(cat3[i].PETRO50FLUX,format='(E8.2)')
      flux90 = string(cat3[i].PETRO90FLUX,format='(E8.2)')
      reff = string(cat2[i].petror90_g,format='(I2)')
      itime = string(cat3[i].petro90itime,format = '(I4)')
      itimetot = string(total_itime(i),format='(I5)')
      if i ge 9 and i lt 99 then nos = string(i+1,format='(I02)')
      if i lt 9 then nos= string(i+1,format='(I01)')
      printf, lun2, nos, ras, decs,gmag, mass,flux50,flux90,reff,itime,itimetot, format='(A2,2(" & ",A11),2(" & ",A5),2(" & ",A10)," & ",A3," & ",A4," & ",A5," \\")'
   
   endfor
    printf,lun2,'\multicolumn{3}total integration time & & & &'+string(total(cat3.petro90itime)/3600.,format='(D4.1)')+'hr &'+string(all_itime,format='(D4.1)')+'hr'
    close,lun2
    free_lun, lun2
    stop
end
