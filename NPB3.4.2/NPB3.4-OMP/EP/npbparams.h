! CLASS = D
!  
!  
!  This file is generated automatically by the setparams utility.
!  It sets the number of processors and the class of the NPB
!  in this directory. Do not modify it by hand.
!  
        character class
        parameter (class ='D')
        integer m
        parameter (m=36)
        logical  convertdouble
        parameter (convertdouble = .false.)
        character compiletime*11
        parameter (compiletime='10 Sep 2023')
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
