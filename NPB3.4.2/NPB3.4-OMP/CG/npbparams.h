! CLASS = C
!  
!  
!  This file is generated automatically by the setparams utility.
!  It sets the number of processors and the class of the NPB
!  in this directory. Do not modify it by hand.
!  
        integer            na, nonzer, niter
        double precision   shift, rcond
        parameter(  na=150000, &
     &              nonzer=15, &
     &              niter=75, &
     &              shift=110., &
     &              rcond=1.0d-1 )
        integer, parameter :: kz=4
        logical  convertdouble
        parameter (convertdouble = .false.)
        character compiletime*11
        parameter (compiletime='09 Sep 2023')
        character npbversion*5
        parameter (npbversion='3.4.2')
        character cs1*8
        parameter (cs1='gfortran')
        character cs2*5
        parameter (cs2='$(FC)')
        character cs3*46
        parameter (cs3='-L/home/ge35vin2/spack/opt/spack/linux-sles...')
        character cs4*46
        parameter (cs4='-I/home/ge35vin2/spack/opt/spack/linux-sles...')
        character cs5*12
        parameter (cs5='-O3 -fopenmp')
        character cs6*9
        parameter (cs6='$(FFLAGS)')
        character cs7*6
        parameter (cs7='randi8')
