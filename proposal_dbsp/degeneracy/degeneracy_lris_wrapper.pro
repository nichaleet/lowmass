pro degeneracy_lris_wrapper
;   if file_test('grid_spec_lris.fits') eq 0 then degeneracy_lris,sn=25,/redo
;   for i=26,39,1 do degeneracy_lris,sn=i
   if file_test('grid_spec_lris_fine.fits') eq 0 then degeneracy_lris,sn=40,/redo,/fine
   for i=20,52 do degeneracy_lris,sn=i,/fine
end
