      PROGRAM ISING_3D
      IMPLICIT NONE
      INTEGER*8:: i,j,L,p,a,b,c,d,mm,nn,qq,N,time,niter,k,count,f,g
      REAL*8:: r,q,E,M,mag,Ei,Ef,dE,u,h
      REAL*8:: T,J_ising=1.0
      REAL*8:: av_m,av_e,cv,chi,av_e2,av_m2,av_m_n,av_e_n
      integer,dimension(:,:,:),allocatable:: spin
      integer:: seed,T_temp,n_stat,n_equil
      character(len=30):: char_a,char_b !char_b stores the name dump_pos
      seed=44859
      char_b="store config"
      print*, "enter the number of lattice points in 1D:"
      read*, L
      print*, "enter the number of iterations:"
      read*, niter
      ALLOCATE(SPIN(L,L,L))
      E=0.0d0 !instantaneous energy of the lattice
      M=0.0d0 !instantaneous magnetic moment
      N=L*L*L !total number of spins
      
      n_equil=10000 ! equilibriation time
      n_stat=10 ! collect statistical data after every n_stat steps
      CALL random_seed
                             
      !-------------------------
      !INITIALIZE THE LATTICE
      !-------------------------
 
      OPEN(21,file="ising_plus.dat")
      OPEN(22,file='ising_minus.dat')
      p=0
      do i=1,L
      do j=1,L
      do k=1,L
      call random_number(r)

      !spin(j,i)=1  ! deterministic initial condition
                   ! All spins points up or down

      if(r.lt.0.5)then
      spin(i,j,k)=-1
      !write(22,*)float(i),float(j)
      else            ! assign spins randomly
      spin(i,j,k)=1
      !write(21,*)float(i),float(j)
      endif
      ! Writing the initial configuration
      !write(21,*)float(i),float(j),float(spin(j,i)) 
      enddo
      enddo
      enddo
      close(21)
      close(22)
      !CALCULATE INITIAL MAGNETIZATION AND ENERGY

      do i=1,L
      do j=1,L
      do k=1,L
      a=i+1
      b=i-1
      c=j+1
      d=j-1
      g=k+1
      f=k-1
      if(i.eq.L)a=1 !PERIODIC BOUNDARY CONDITION
      if(i.eq.1)b=L
      if(j.eq.1)d=L
      if(j.eq.L)c=1
      if(k.eq.L)g=1
      if(k.eq.1)f=L
      
       E=E-J_ising*dfloat((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)  &
         +spin(i,c,k)+spin(i,d,k)+spin(i,j,g)+spin(i,j,f)))
      M=M+spin(i,j,k)

      enddo
      enddo
      enddo
      mag=M/(dfloat(N))  ! total number of spins: L*L
      E=E*0.5         ! to avoid double counting
      print*, "Initial energy E, energy per spin:",E,E/DFLOAT(N)
      print*, "Initial marnetization, mag. per spin:",M,mag

      ! INITIALIZATION DONE
      ! PROCEED TO REACH EQUILIBRIUM

      open(unit=10,file="data.dat")    
      do T_temp=470,380,-2             ! initiation of the temperature loop
      T=dfloat(T_temp)/100.0d0        ! fix T  
      
      av_m=0.0d0
      av_e=0.0d0
      av_m_n=0.0d0
      av_e_n=0.0d0
      av_m2=0.0d0
      av_e2=0.0d0
      
      do count=1,niter
      do mm=1,L
      do nn=1,L
      do qq=1,L
      
       call random_number(r); i=int(r*float(L))+1 !random lattice site
       call random_number(r); j=int(r*float(L))+1
       call random_number(r); k=int(r*float(L))+1

       !********************************************************
       ! ASSIGNMENT OF NEIGHBOURING SITES FOR EACH LATTICE POINT
       !********************************************************

       a=i+1
       b=i-1
       c=j+1
       d=j-1
       g=k+1
       f=k-1

      !******************************************
      ! IMPOSITION OF PERIODIC BOUNDARY CONDITION
      !******************************************

      if(i.eq.L)a=1 !PERIODIC BOUNDARY CONDITION
      if(i.eq.1)b=L
      if(j.eq.1)d=L
      if(j.eq.L)c=1
      if(k.eq.L)g=1
      if(k.eq.1)f=L


      !************************************
      ! PERIODIC BOUNDARY CONDITION IMPOSED
      !************************************

      !****************************************************       
      ! CHECK THE SIGN OF ENERGY SHIFT BECAUSE OF SPIN FLIP
      !****************************************************
      
       Ei=E-J_ising*dfloat((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)  &
         +spin(i,c,k)+spin(i,d,k)+spin(i,j,g)+spin(i,j,f)))

      spin(i,j,k)=-spin(i,j,k)  ! trial flip
      Ef=E-J_ising*dfloat((spin(i,j,k))*(spin(a,j,k)+spin(b,j,k)  &
         +spin(i,c,k)+spin(i,d,k)+spin(i,j,g)+spin(i,j,f)))
      dE=Ef-Ei


      !===================================
      ! IMPOSITION OF METROPOLIS ALGORITHM
      !===================================


      if(dE.le.0.0d0)then
        E=E+dE
        M=M+2.0d0*dfloat(spin(i,j,k))    
      else
        u=exp(-dE/T)
        call random_number(h)
      if(h.lt.u)then
        E=E+dE
        M=M+2.0d0*dfloat(spin(i,j,k))  ! parameters updated
      else
        spin(i,j,k)=-spin(i,j,k) ! flip not accepted; parameters not updated
      endif
      endif
      enddo
      enddo
      enddo

      !======================================================
      ! AFTER REACHING EQUILIBRIUM, COLLECT STATISTICAL DATA
      !======================================================
      
      if(count.gt.n_equil)then
      !if(mod(k,n_stat).eq.0)then
      mag=M/dfloat(N)
      av_m=av_m+mag
      av_e=av_e+E/dfloat(N)
      av_m_n=av_m_n+M
      av_e_n=av_e_n+E
      av_m2=av_m2+(M*M)    
      av_e2=av_e2+(E*E)
      !endif
      endif
      enddo

      av_m=av_m/dfloat(niter-n_equil)
      av_e=av_e/dfloat(niter-n_equil)
      av_e2=av_e2/dfloat(niter-n_equil)
      av_e_n=av_e_n/dfloat(niter-n_equil)
      av_m2=av_m2/dfloat(niter-n_equil)
      av_m_n=av_m_n/dfloat(niter-n_equil)
      
      cv=(av_e2-av_e_n*av_e_n)/T*T 
      chi=(av_m2-av_m_n*av_m_n)/T 
      
      write(10,*)T,av_m,av_e,cv,chi
      enddo 

      close(10)
      deallocate(spin)
      END PROGRAM ISING_3D
