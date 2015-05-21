       subroutine input
        
       use technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'


       integer :: i,il,ij,lll,ak             !single particle input
       integer :: a,b,cc,c,ff,ss              !phonons input
       integer :: ind,par,jj,inum,dimea
       integer :: int1,int2,int3,int4
       integer, allocatable :: ll(:),part(:,:)
       double precision :: en,val
       double precision :: e,qe,v     !single particle input

!**********************************************************************
!******************input file******************************************
!**********************************************************************

      open(1,file='input.dat',status='old',form='formatted')   
   
      read(1,*)AN,AZ
      read(1,*)ihnmin,ihnmax
      read(1,*)ihpmin,ihpmax
      read(1,*)ipnmin,ipnmax
      read(1,*)ippmin,ippmax
      read(1,*)iparmin,iparmax
      read(1,*)jjmin,jjmax
      read(1,*) iparmin2,iparmax2
      read(1,*) jjmin2,jjmax2
      read(1,*)iv
      close(1)
       write(*,*) 'hole proton coup to phonon problem'
       write(*,*) 'A=',AN,'Z=',AZ
       write(*,*) 'ho=',hbarom,'MeV'
       write(*,*) 'TDA phonons from j=',jjmin,'to',jjmax
       write(*,*) 'phonon parity from',iparmin,iparmax
       write(*,*) 'Proton space',ihpmin,ihpmax,ippmin,ippmax
       write(*,*) 'Neutron space',ihnmin,ihnmax,ipnmin,ipnmax
       write(*,*) '2 phonon parity ',iparmin2,iparmax2
       write(*,*) '2 phonon J from',jjmin2,jjmax
       write(*,*) 'iv=',iv
!**********************************************************************
!********************single particle levels****************************
!**********************************************************************    

       open(12,file='../TDA/HF_n.out',status='old',form='formatted')
        read(12,*)

       do while (.not.eof(12))
        read(12,*)i,il,ij,e,qe,v
       enddo
        imax=i
       allocate(levn(imax))
       close(12)
         open(12,file='../TDA/HF_n.out',status='old',form='formatted')
          read(12,*)
          do while (.not.eof(12))
         read(12,*)i,il,ij,e,qe,v
          levn(i)%index=i
          levn(i)%l=il
          levn(i)%j2=ij
          levn(i)%ei=e
        enddo

       close(12)

       id=imax

       open(12,file='../TDA/HF_p.out',status='old',form='formatted')
        read(12,*)

        do while (.not.eof(12))
        read(12,*)i,il,ij,e,qe,v
        enddo
       allocate(levp(i))
        close(12)

         open(12,file='../TDA/HF_p.out',status='old',form='formatted')
         read(12,*)

        do while (.not.eof(12))
        read(12,*)i,il,ij,e,qe,v
        levp(i)%index=i
        levp(i)%l=il
        levp(i)%j2=ij
        levp(i)%ei=e
        enddo
        if (i.gt.imax) imax=i
        close(12)
       jmax=0

        do i=1,id
         if(levp(i)%j2.gt.jmax) jmax=levp(i)%j2
        enddo
         do i=1,id
          if(levn(i)%j2.gt.jmax) jmax=levn(i)%j2
        enddo

       Jtot=levp(iv)%j2
       ipartot=(-1)**levp(iv)%l

       write(*,*) 'J e parità dello stato da generare'
       write(*,*)  Jtot,ipartot
!*************************************************************************
!********************phonons**********************************************
!*************************************************************************

      write(*,*) jjmin,jjmax,jmax

      open(1,file='../TDA/phon_inf.out',status='old',form='formatted')
       num_ph=0 
       read(1,*) c
       allocate(dime(c))
       numb=c
       do is=1,c
       read(1,*) a,b,cc
       dime(is)%parity=a
       dime(is)%jph=b
       dime(is)%dimen=cc
       num_ph=num_ph+cc 
        end do

       close(1)
         write(*,*)'num_ph=',num_ph
       allocate(ph(num_ph))
!*******************************************************************************
        open(3,file='../TDA/list.out',status='old',form='unformatted')
          inum=0
          do while(.not.eof(3))
          inum=inum+1

         read(3) ind,par,jj,en

         ph(inum)%par=par
         ph(inum)%Jph=jj
         ph(inum)%lam=ind
         ph(inum)%el=en 
          end do

         close(3)
!***********************************************************************************
!***************************************2p-ph basis*********************************
!***********************************************************************************
       open(4,file='infobase.dat',status='old',form='formatted')
      
       dimensi=0
       do while(.not.eof(4))
        read(4,*) as,ae,ar
         dimensi=dimensi+1
        end do
        close(4)
         dimbas=dimensi
         write(*,*) 'dimension of coup. basis=',dimbas
       open(5,file='infobase.dat',status='old',form='formatted')
          allocate(bas(dimbas))
          do while(.not.eof(5))
         read(5,*) ii,ir,ilam
         bas(ii)%ir=ir
         bas(ii)%ilam=ilam
          end do

         close(5)
!************************************************************************************
!**********************************************************************
      open(12,file='2ph/2f_x.dat',status='old',form='unformatted')
      allocate(dime2ph(-1:1,0:jjmax2))
      ss=0
      dime2ph=0
      do while(.not.eof(12))
       read(12) int1,int2,int3,int4
        dime2ph(int1,int2)=int3
         do j=1,int3
         read(12)(ila,ilap,tt,i=1,int4)
          end do
         ss=ss+int3
         end do
       close(12)
       write(*,*) ss,jjmax2
       allocate(ll(0:100))
       allocate(part(-1:1,0:jjmax2))
       ll=0
       part=0
       tt=0
       do ii=-1,1
        do jj=0,jjmax2
          if(ii.ne.0) then
          tt=tt+1
         ll(0)=0
         ll(tt)=ll(tt-1)+dime2ph(ii,jj)
         part(ii,jj)=ll(tt-1)
         end if
        end do
       end do
!************************************************************************
      open(1,file='2ph/2f_states.dat',status='old',form='unformatted')
         dimea=0
         do while(.not.eof(1))
          read(1) int1,int2,int3,en
          dimea=dimea+1
          end do 
        close(1)
       num_2ph=dimea
          ll=0
      open(1,file='2ph/2f_states.dat',status='old',form='unformatted')
          allocate(ph2(dimea))
         do while(.not.eof(1))
          read(1) int1,int2,int3,en
    
          ph2(int1)%par=int2
          ph2(int1)%J=int3
          ph2(int1)%el=en
          ph2(int1)%alfa=int1-part(int2,int3)
          end do
        close(1)
        write(*,*)'num_2ph=',num_2ph
!************************************************************************************
       open(6,file='dmat.dat',status='old',form='formatted')
        allocate(dmat(dimbas,dimbas))
       dmat=0.d0
       do while(.not.eof(6))
        read(6,*) row,col,dd
        dmat(row,col)=dd
   
       end do
       close(6)
!************************************************************************************* 
      open(7,file='d1matr.dat',status='old',form='formatted')
      read(7,*) ak,no
      write(*,*) ak,no
      allocate(d1(no,ak))
      do while(.not.eof(7))
      read(7,*) i,j,d1(i,j)
      end do
  
      close(7)
!**************************************************************************************
      allocate(mxt(dimbas),mxt1(dimbas))  
      open(8,file='mixt.dat',status='old',form='formatted')

      do while(.not.eof(8))
      read(8,*)i,mm
       mxt(i)=mm
       mxt1(mm)=i
      enddo

      close(8)

!**************************************************************

      open(9,file='eigenDAD.dat',status='old',form='formatted')
      allocate(eigen(no))
      do while(.not.eof(9))
      read(9,*)i,w
       eigen(i)=w
      enddo

      close(9)
!***************************************************************
       Allocate(eigv(no,no))
       eigv=0.d0
       open(44,file='eigvDAD.dat',status='unknown',form='formatted')
        do while(.not.eof(44))
       read(44,*) ii,jj,val
        eigv(ii,jj)=val
        end do
        close(44)

!***************************************************************


      return
        end subroutine
