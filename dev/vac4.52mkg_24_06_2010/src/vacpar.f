*##############################################################################
* include vacpar - mhd

* For MHD: density,momentum,energy,magnetic_field + adiabatic_index,resistivity

*     VACPHYS module name
      CHARACTER*3  typephys
      PARAMETER(  typephys='mhd')
*     Equation parameter names
      CHARACTER*9  eqparname
      PARAMETER(  eqparname='gamma eta')

*     flow variables
      INTEGER  rho_,m0_,m1_,m2_,e_
      PARAMETER(  rho_=1,m0_=rho_,m1_=m0_+1,m2_=m0_+2,e_=m2_+1)
*     flow variables
      INTEGER  ee_, b0_,b1_,b2_
      PARAMETER(  ee_=e_, b0_=e_,b1_=b0_+1,b2_=b0_+2)
*     No. flow variables
      INTEGER  nw
      PARAMETER(  nw=b2_)

*     Primitive variables
      INTEGER  v0_, v1_,v2_, p_, pp_
      PARAMETER(  v0_=m0_, v1_=m1_,v2_=m2_, p_=e_, pp_=ee_)

*     Characteristic
      INTEGER  fastRW_,fastLW_,slowRW_,slowLW_
      PARAMETER(  fastRW_=1,fastLW_=2,slowRW_=3,slowLW_=4)
*     waves
      INTEGER  entroW_,diverW_,alfvRW_,alfvLW_
      PARAMETER(  entroW_=5,diverW_=6,alfvRW_=7,alfvLW_=8)
*     Potential extrema
      INTEGER  extremeRW_,extremeLW_
      PARAMETER(  extremeRW_=fastRW_,extremeLW_=fastLW_)

*     Polar var. names
      INTEGER  mr_,mphi_,mz_
      PARAMETER(  mr_=m0_+r_,mphi_=m0_+phi_,mz_=m0_+z_)
      INTEGER  br_,bphi_,bz_
      PARAMETER(  br_=b0_+r_,bphi_=b0_+phi_,bz_=b0_+z_)

*     No. vector vars
      INTEGER  nvector
      PARAMETER(  nvector=2)

*     equation params
      INTEGER  gamma_,eta_,neqpar
      PARAMETER(  gamma_=1,eta_=2,neqpar=2)
      INTEGER  divbcoeff_,divbconst_,divbbound_,nprocpar
      PARAMETER(  divbcoeff_=1,divbconst_=2,divbbound_=3,nprocpar=3)
*                                                               processing params

* end include vacpar - mhd
*##############################################################################
