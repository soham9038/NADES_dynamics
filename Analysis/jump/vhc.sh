## Script to calculate the van Hove correlation ##


for i in {1..66..1} # We have 66 glucose molecules in the system
do
cat > vhc.f90<<EOF      
      Program TAHCF
      Implicit None
 
      Integer, Parameter     :: DP = KIND(1.d0)
      Integer, Parameter:: Nstep = 9996121, k = 4517  ! Nstep = k * number of segments
      Integer ::  i, j
      character*40 :: input, output
      Real*4 :: dRx, dRy, dRz, dt, dR!, t
      Real*4 :: Rx(Nstep), Ry(Nstep), Rz(Nstep)

3     format(F20.16)     ! 16 Digits in numbers 

      !Open (unit = 12, file = 'name.dat', status = 'OLD') 
      !Read(12,*) input
      !Read(12,*) output
      Open (Unit = 13, File = 'com_unfold_${i}.xyz', Status = 'OLD') ! Input file
      Open (Unit = 15, File = 'VHC_${i}.xyz', Status = 'UNKNOWN') ! Output file


      do i = 1, Nstep

      Read(13,*) Rx(i), Ry(i), Rz(i)

      end do
     
      do j = 1, Nstep-k
       
         
           dRx = Rx(j+k) - Rx(j)
	   dRy = Ry(j+k) - Ry(j)
           dRz = Rz(j+k) - Rz(j)
           dR = sqrt(dRx*dRx + dRy*dRy + dRz*dRz) 

           
      write(15,3) dR
      
      enddo
      
      end program TAHCF
EOF
gfortran vhc.f90
./a.out
done

cat VHC_*.xyz > combined.xyz
