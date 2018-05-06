pro get_kcorrectmass,str

   cgals = n_elements(str)
   ;magnitudes are in AB maggie sdss system 
   color = dblarr(5,cgals);u,g,r,i,z
   color[0,*] = str.PETROMAG_U
   color[1,*] = str.PETROMAG_G
   color[2,*] = str.PETROMAG_R
   color[3,*] = str.PETROMAG_I
   color[4,*] = str.PETROMAG_Z

   color_err = dblarr(5,cgals) 
   color_err[0,*] = str.PETROMAGERR_U
   color_err[1,*] = str.PETROMAGERR_G
   color_err[2,*] = str.PETROMAGERR_R
   color_err[3,*] = str.PETROMAGERR_I
   color_err[4,*] = str.PETROMAGERR_Z
   color_ivar = 1/color_err^2
   filterlist=['sdss_u0.par', 'sdss_g0.par','sdss_r0.par', 'sdss_i0.par','sdss_z0.par']
   zlist = str.z

   kcorrect,color,color_ivar,zlist,kcorrect1,chi2=chi2,filterlist=filterlist,/sdssfix,mass=mass,b300=b300
   mass = mass/0.49 ;fix for the hubble constant h=0.7
   mass = alog10(mass)

   str.KCORRECTMASS = mass

   color = dblarr(5,cgals);u,g,r,i,z
   color[0,*] = str.DERED_U
   color[1,*] = str.DERED_G
   color[2,*] = str.DERED_R
   color[3,*] = str.DERED_I
   color[4,*] = str.DERED_Z

   color_ivar = 1/color_err^2
   filterlist=['sdss_u0.par', 'sdss_g0.par','sdss_r0.par', 'sdss_i0.par','sdss_z0.par']
   zlist = str.z

   kcorrect,color,color_ivar,zlist,kcorrect1,chi2=chi2,filterlist=filterlist,/sdssfix,mass=mass,b300=b300
   mass = mass/0.49 ;fix for the hubble constant h=0.7
   mass = alog10(mass)

   str.KCORRECTMASS_DERED = mass
end
