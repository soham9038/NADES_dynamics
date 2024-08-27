## Script to calculte RMSD and LAMBDA ##

#com file for individual molecule
for m in {1..111..1}
do
	echo "$m" > input.txt
	gmx traj -f md.xtc -s md.tpr -b 49900 -e 50000 -n sol.ndx -com -ox com_"$m".xvg < input.txt
	paste com_"$m".xvg | tail -n10001 > test.txt
	paste test.txt | awk '{print $2, $3, $4}' > com_"$m".xyz
	rm -rf com_"$m".xvg
done

#unfolding considering singleparticle code
for i in {1..111..1}
do
cat > UNFOLD.f90<<EOF
      Program unfold
      Implicit None
 
      Integer, Parameter     :: DP = KIND(1.d0)
      Integer, Parameter:: Nstep = 10001
      Integer ::  i
      Real*4 :: dRx, dRy, dRz, dR
      Real*4 :: Rx(Nstep), Ry(Nstep), Rz(Nstep) !t(Nstep)
      character :: t(Nstep)
      Open (unit = 13, file = 'com_${i}.xyz', status = 'OLD') 
      Open (Unit = 15, File = 'com_unfold_${i}.xyz', Status = 'UNKNOWN')


      do i = 1, Nstep

      Read(13,*) Rx(i), Ry(i), Rz(i)

      end do

      do i = 2, Nstep
           dRx = Rx(i) - Rx(i-1)
           dRy = Ry(i) - Ry(i-1)
           dRz = Rz(i) - Rz(i-1)
           dRx = dRx - ANINT(dRx/4.23267)*4.23267
           dRy = dRy - ANINT(dRy/4.23267)*4.23267
           dRz = dRz - ANINT(dRz/4.23267)*4.23267
           Rx(i) = Rx(i-1) + dRx
           Ry(i) = Ry(i-1) + dRy
           Rz(i) = Rz(i-1) + dRz

       write(15, *) Rx(i), Ry(i), Rz(i)

      end do

      end program unfold
EOF
gfortran UNFOLD.f90
./a.out
done

#COM, LAMBDA, Rg calculation
for i in {1..111..1}
do
cat > rg.f90<<EOF      
     Program COM_Radius_Gyration
     Implicit None
 
      Integer, Parameter     :: DP = KIND(1.d0)
      Real*4, Allocatable, Dimension(:) :: P(:), Q(:), M(:), Rgx(:), Rgy(:), Rgz(:), RG(:)
      Integer, Parameter:: Nstep = 9999248, Np = 1168 
      Integer ::  i, j, k, l
      Real*4 :: Rx(Nstep), Ry(Nstep), Rz(Nstep), t(Nstep)
      Allocate(P(Nstep), Q(Nstep), M(Nstep), Rgx(Nstep), Rgy(Nstep), Rgz(Nstep), RG(Nstep))

      Open (Unit = 12, File = 'com_unfold_{i}.xyz', Status = 'OLD')
      Open (Unit = 15, File = 'COM_{i}.dat', Status = 'UNKNOWN')
      Open (Unit = 18, File = 'RG_{i}.dat', Status = 'UNKNOWN') ! RMSD file
      Open (Unit = 20, File = 'Lambda_{i}.dat', Status = 'UNKNOWN') ! Lambda file
      
      
      
      do i = 1, Nstep

      Read(12, *) Rx(i), Ry(i), Rz(i)

      end do
      
      do j = 1, Nstep, Np
      
         P(j) = 0.0d0
         Q(j) = 0.0d0
         M(j) = 0.0d0

        do i= 1, Np
       

          P(j)= (P(j) + Rx(i+j-1))  

          Q(j)= (Q(j) + Ry(i+j-1))  

          M(j)= (M(j) + Rz(i+j-1)) 
     
        end do
         
          write(15,*) P(j)/Np, Q(j)/Np, M(j)/Np
          
      end do
      
      
      do k = 1, Nstep, Np
      
         Rgx(k) = 0.0d0
         Rgy(k) = 0.0d0
         Rgz(k) = 0.0d0
         RG(k)  = 0.0d0
           
           do l = 1, Np
           
              Rgx(k) = Rgx(k) + (Rx(l+k-1) - (P(k)/Np))**2
              
              Rgy(k) = Rgy(k) + (Ry(l+k-1) - (Q(k)/Np))**2
              
              Rgz(k) = Rgz(k) + (Rz(l+k-1) - (M(k)/Np))**2
              
           end do
           
          RG(k) = sqrt((Rgx(k) + Rgy(k) + Rgz(k))/Np)
          
          write(18,*)  RG(k)
          write(20,*) 2*RG(k)
      
      end do 
      
      Deallocate(P, Q, M, Rgx, Rgy, Rgz, RG)
      CLOSE(12)
      CLOSE(15)
      CLOSE(18)
         
      end program COM_Radius_Gyration

EOF
gfortran rg.f90
./a.out
done

