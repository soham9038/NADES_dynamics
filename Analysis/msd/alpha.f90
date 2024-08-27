!! Fortran code to calculate the time dependence exponent of MSD curve !!


      Program Alpha
      Implicit None
 
      Integer, Parameter     :: DP = KIND(1.d0)
      Real*4, Allocatable, Dimension(:) :: dy(:), dx(:), A(:)
      Integer, Parameter:: Nstep = 10000001, m = 50 ! Nstep = Total number of stes involved, m = number steps skipped during calculation
      Integer ::  i
      Real*8 :: Rx(Nstep), t(Nstep)
      Allocate(dy(Nstep), dx(Nstep), A(Nstep))


       Open (Unit = 13, File = 'msd_glc.dat', Status = 'OLD') !input file name
       Open (Unit = 15, File = 'alpha_glc_50p.dat', Status = 'UNKNOWN') ! output file name


      do i= 1, Nstep
        
        Read(13,*) t(i), Rx(i) 

      end do

      do i= 1, Nstep-m

         dy(i) = (log(Rx(i+m)))-(log(Rx(i)))
         dx(i) = (log(t(i+m)))-(log(t(i)))

         A(i) = dy(i)/dx(i)

        write(15,*) t(i), A(i)
      end do

      end program Alpha
    
