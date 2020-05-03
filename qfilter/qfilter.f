c      This program, QFILTER modifies classical trajectories according 
c      to the GSTA method:
c      Dénes Berta, Dávid Ferenc, Imre Bakó and Ádám Madarász 
c      "Nuclear Quantum Effects from the Analysis of Smoothed Trajectories: 
c      Pilot Study for Water"
c      https://doi.org/10.1021/acs.jctc.9b00703
c      email: madarasz.adam@ttk.hu
c     
c      QFILTER works interactvely, just start the binary file, then
c      1. you need to give the name of the xyz trajectory file
c      2. the program asks for the timestep between two snapshots
c      3. the temperature of the simulation must be entered
c      

      program qfilter

c These parameters determines the numerical integrals:
      DOUBLE PRECISION ftstep,kercut,precint

      parameter (ftstep=1.0d-4)    ! step size in the Fourier transform
      parameter (kercut=1.0d-6)    ! cutoff value of the kernel function
      parameter (precint=1.0d-11)  ! precision of the integral

      DOUBLE PRECISION avogadro,boltzmann,planck

      parameter (avogadro=6.02214129d+23)
      parameter (boltzmann=0.831446215d0)
      parameter (planck=6.62606957d-34)

      integer ndec,ablak,nblock,natoms
      integer i,j,k,m,ipos,tav,midpos
      DOUBLE PRECISION kernel(0:100000),PI,temp
      DOUBLE PRECISION dnu,sum,koz,timestep,deltat,v,seg,ana
      DOUBLE PRECISION, allocatable :: xtot(:,:),ytot(:,:),ztot(:,:)
      DOUBLE PRECISION, allocatable :: cx(:),cy(:),cz(:)
      character*120 xyzfile,outxyzfile
      character(75), allocatable :: title(:)  
      character(5), allocatable :: nam(:)  

      PI=4.D0*DATAN(1.D0)

c
c     get the base name of user specified input structures
c
      write (*,10)
  10    format (/,' Enter the name of the xyz file: ',$)
      read (*,20)  xyzfile
   20    format (a120)

c
c     get the number of frames
c

      OPEN (50, file = xyzfile)
      read(50,*) natoms
      i = 1
      DO
       READ (50,*, END=11)
       i = i + 1
      END DO
   11 CLOSE (50)

      nblock=i/(natoms+2)

      write(*,*) 'Number of atoms: ',natoms
      write(*,*) 'Number of frames: ',nblock

      write (*,80)
   80 format (/,' Timestep between frames',
     &              ' in picoseconds:  ',$)
   90 format (f20.0)
      read (*,90) timestep

      write (*,81)
   81 format (/,' Temperature for filtration',
     &             ' in Kelvin:  ',$)
   91    format (f20.0)
         read (*,91) temp

      write(*,*) 'Calculation of the kernel function:'

      deltat=2*PI/(avogadro*planck/boltzmann/temp/timestep*1E11)

      j=0
      sum=1.0

      do while (abs(sum) .gt. kercut)

         koz=0.5

         i=0
         sum=0.0
         dnu=ftstep/(dble(j+1))

         do while (koz .gt. precint)
            sum=sum+koz*cos(v*dble(j)*deltat)
            i=i+1
            v=dnu*dble(i)
            koz=dsqrt(v/2.0)*(1.0/dsqrt(tanh(v/2.0))-1.0)
         enddo
      
         seg=sum*dnu/PI*deltat
    
         ana=0.5/dsqrt(PI*(dble(j)+0.5)*deltat)
      ana=ana-sign(0.5,float(j)-0.5)/dsqrt(PI*(abs(dble(j)-0.5))*deltat)
         sum=seg+ana

         write(*,*) sum
         kernel(j)=sum
         j=j+1
      enddo

      write(*,*) 'Starting filtration'
c     End of the calculation of the kernel function

      ndec=j-1
      ablak=2*ndec+1

c
c     perform dynamic allocation of some local arrays
c
      allocate (title(nblock))
      allocate (nam(natoms))
      allocate (cx(natoms))
      allocate (cy(natoms))
      allocate (cz(natoms))

      allocate (xtot(natoms,ablak))
      allocate (ytot(natoms,ablak))
      allocate (ztot(natoms,ablak))

      open (unit=50,file=xyzfile,status='old',action='read')  

      outxyzfile='filt_'//xyzfile

      open (unit=60,file=outxyzfile,status='unknown',action='write')  

c
c     cycle over all pairs of snapshot frame blocks
c
      do i = 1, nblock

         read(50,*) natoms
         read(50,15) title(i)
   15    format (a100)

         do j=1,natoms
            read(50,*) nam(j),cx(j),cy(j),cz(j)  
         end do  
 
         ipos=mod(i-1,ablak)+1
         midpos=mod(i-1-ndec,ablak)+1
         do j = 1, natoms
            xtot(j,ipos)=cx(j)
            ytot(j,ipos)=cy(j)
            ztot(j,ipos)=cz(j)
         end do

         if (i .ge. ablak) then

            do m=1, natoms
               cx(m)=0.0d0
               cy(m)=0.0d0
               cz(m)=0.0d0
            end do
               
            do k=1,ablak
               tav=min(abs(midpos-k),ablak+k-midpos,ablak+midpos-k)
               do m=1, natoms
                  cx(m)=cx(m)+xtot(m,k)*kernel(tav)
                  cy(m)=cy(m)+ytot(m,k)*kernel(tav)
                  cz(m)=cz(m)+ztot(m,k)*kernel(tav)
               end do
            end do

c
c     write output
c
      
            write(60,*) natoms
            write(60,*) title(i-ndec)

            do k=1,natoms
               write(60,*) nam(k),cx(k),cy(k),cz(k)  
            end do  
			
         end if
         
      end do

c
c     perform deallocation of some local arrays
c
      deallocate (title)
      deallocate (nam)
      deallocate (cx)
      deallocate (cy)
      deallocate (cz)

      deallocate (xtot)
      deallocate (ytot)
      deallocate (ztot)

      close(50)
      close(60)

      write(*,*) 'Number of filtered frames: ', nblock-ablak+1
      write(*,*) 'Name of the outputfile: ', outxyzfile

      end program qfilter
