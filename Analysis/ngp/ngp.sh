## Calculation of non-gaussian parameter ##

#centre-of-mass calculation for individual molecules

for m in {1..133..1} # We have 133 urea molecules in the system
do
	echo "$m" > input.txt
	gmx traj -f md.xtc -s md.tpr -b 90000 -e 100000 -n ../ure.ndx -com -ox com_"$m".xvg < input.txt # 90 ns to 100 ns is chosen as representative block. Needs to be averaged over whole trajectory
	paste com_"$m".xvg | tail -n1000001 > test.txt
	paste test.txt | awk '{print $2, $3, $4}' > com_"$m".xyz
	rm -rf com_"$m".xvg
done


# unfolding individual molecular trajectory #

# Nstep: Number of steps
# input file: com_${i}.xyz
# output file: com_unfold_${i}.xyz
# box edge length: 2.97605 nm

for i in {1..133..1}
do
cat > UNFOLD.f90<<EOF
      Program unfold
      Implicit None
 
      Integer, Parameter     :: DP = KIND(1.d0)
      Integer, Parameter:: Nstep = 1000001 
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
           dRx = dRx - ANINT(dRx/2.97605)*2.97605
           dRy = dRy - ANINT(dRy/2.97605)*2.97605
           dRz = dRz - ANINT(dRz/2.97605)*2.97605
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

# ngp calculation #

# 133 urea molecules
# Ntau = maximum lagtime
# input file: com_unfold_${i}.xyz
# output file: NGP_ure_${i}.xyz

for i in {1..133..1}
do
cat > ngp.f90<<EOF      
      Program NGP
      Implicit None
 
      Integer, Parameter     :: DP = KIND(1.d0)
      Real*4, Allocatable, Dimension(:) ::  tau(:), N(:), De(:), Nu(:)
      Integer, Parameter:: Nstep = 1000000
      Integer, Parameter:: Ntau = 100001
      Integer ::  i, j
      character*40 :: input, output
      Real*4 :: dRx, dRy, dRz
      Real*4 :: Rx(Nstep), Ry(Nstep), Rz(Nstep)!, Tot
      Allocate(tau(Ntau), N(Ntau), De(Ntau), Nu(Ntau) )

      
#      !Open (unit = 12, file = 'name.dat', status = 'OLD') 
#      !Read(12,*) input
#      !Read(12,*) output
      Open (Unit = 13, File = 'com_unfold_${i}.xyz', Status = 'OLD')
      Open (Unit = 15, File = 'NGP_ure_${i}.xyz', Status = 'UNKNOWN')

      do i = 1, Nstep

      Read(13,*) Rx(i), Ry(i), Rz(i)

      end do

      do j = 1, Ntau
     
         tau(j) = (j-1)

         Nu = 0.0d0
         De = 0.0d0
      do i = 1, Nstep-j+1
       
            dRx = Rx(i+j-1) - Rx(i)
            dRy = Ry(i+j-1) - Ry(i)
            dRz = Rz(i+j-1) - Rz(i)
            De(j) = De(j) + (dRx*dRx + dRy*dRy + dRz*dRz)
           
            Nu(j)= Nu(j) + ((dRx*dRx + dRy*dRy + dRz*dRz)*(dRx*dRx + dRy*dRy + dRz*dRz))
            

      end do
      
       De(j) = 5*((De(j) / (Nstep-j))*(De(j) / (Nstep-j)))
       Nu(j) = 3*(Nu(j) / (Nstep-j))
       N(j) = (Nu(j) / De(j)) - 1

       write(15, *) (tau(j)), N(j)

      end do

      end program NGP
EOF
gfortran ngp.f90
./a.out
done

# Post processing #

awk '{a[FNR]+=$1;b[FNR]+=$2} END{for (i=1;i<=FNR;i++) print a[i], b[i]}' NGP_ure_*.xyz > file.txt
awk '{print ($1)/133, ($2)/133}' file.txt > ngp.txt
