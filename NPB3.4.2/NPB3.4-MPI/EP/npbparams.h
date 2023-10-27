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
        character*11 compiletime
        parameter (compiletime='25 Aug 2023')
        character*5 npbversion
        parameter (npbversion='3.4.2')
        character*6 cs1
        parameter (cs1='mpif90')
        character*8 cs2
        parameter (cs2='$(MPIFC)')
        character*46 cs3
        parameter (cs3='-L/home/ge35vin2/spack/opt/spack/linux-sles...')
        character*46 cs4
        parameter (cs4='-I/home/ge35vin2/spack/opt/spack/linux-sles...')
        character*29 cs5
        parameter (cs5='-O3 -fallow-argument-mismatch')
        character*9 cs6
        parameter (cs6='$(FFLAGS)')
        character*6 cs7
        parameter (cs7='randi8')
