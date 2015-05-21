       subroutine input_files
        
       use technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       integer :: i,il,ij             !single particle input
       integer :: a,b,cc,c               !phonons input
       integer :: ind,par,jj,inum

       double precision :: en
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
      read(1,*)hbarom
      read(1,*)iparmin,iparmax
      read(1,*)jjmin,jjmax
      read(1,*)iv
      close(1)

       if(iv.eq.0) then
       write(*,*) '2j state to generate?'
       read(*,*)  Jtot
       write(*,*) 'parity of the state?'
       read(*,*)  ipartot
       end if

       write(*,*) 'calculation for proton hole!'

       write(*,*) 'A=',AN,'Z=',AZ
       write(*,*) 'ho=',hbarom,'MeV'
       write(*,*) 'TDA phonons from j=',jjmin,'to',jjmax
       write(*,*) 'phonon pèarity from',iparmin,iparmax
       write(*,*) 'Proton space',ihpmin,ihpmax,ippmin,ippmax
       write(*,*) 'Neutron space',ihnmin,ihnmax,ipnmin,ipnmax
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
        levn(i)%ipar=(-1)**levn(i)%l
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
        levp(i)%ipar=(-1)**levp(i)%l
        enddo

       if (i.gt.imax) imax=i

       jmax=0

       do i=1,id
        if(levp(i)%j2.gt.jmax) jmax=levp(i)%j2
       enddo
       do i=1,id
        if(levn(i)%j2.gt.jmax) jmax=levn(i)%j2
       enddo

       !Jtot=levn(iv)%j2
       !ipartot=(-1)**levn(iv)%l

       write(*,*) 'J e parità dello stato da generare'
       write(*,*)  Jtot,ipartot
!*************************************************************************
!********************phonons**********************************************
!*************************************************************************
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
         return
        end subroutine
