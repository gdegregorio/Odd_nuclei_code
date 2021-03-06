       program Hamiltonian
  

       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
       
       double precision, allocatable, save ::ham(:,:),wro(:)
       double precision, allocatable, save :: wr(:),wi(:),vr(:,:)
       double precision, allocatable, save :: dmat1(:,:),term(:)
       double precision :: vec

     
      call input
      call Finteraction
      call Vinteraction

      allocate(dmat1(dimbas,dimbas))
      dmat1=dmat

      call permut(dimbas,dmat,mxt)
      !write(*,*) 'perm dmat'
      call xmatr
      !write(*,*) 'xma'
      call permutv(dimbas,Fcoup,mxt)
      !write(*,*) 'perm Fco'
       allocate(term(no))  
        term=0.d0
        open(897,file='coup_term.dat',status='unknown',form='formatted')

       do j=1,no
        vec=0.d0
 
        do i=1,dimbas
         Jl=ph(bas(i)%ilam)%Jph
         vec=vec+(dble(2*Jl+1))**0.5d0*Fcoup(i)*xm(i,j)
        end do
        term(j)=vec!(1.d0/dsqrt(dble(levn(iv)%j2)+1))*vec
        write(897,*) j,term(j)
      
        end do
         close(897)

        allocate(ham(no+1,no+1))
         ham=0.d0
        ham(1,1)=levn(iv)%ei
        do i=2,no+1
          ham(i,1)=term(i-1)
          ham(1,i)=term(i-1)
          ham(i,i)=eigen(i-1)
        end do
        open(312,file='DAD+coup.dat',status='unknown',form='formatted')
       do i=1,no+1
        do j=1,no+1
         write(312,*) i,j,ham(i,j)
        end do
       end do
        close(312)
       allocate(wr(no+1),wi(no+1),vr(no+1,no+1))
 
       wr=0.d0
       wi=0.d0 
       vr=0.d0
  
      call diagonalization(ham,wr,wi,vr,no+1)
      
      open(334,file='eigencoup.dat',status='unknown',form='formatted')
      open(335,file='ordeigcoup.dat',status='unknown',form='formatted')
       write(334,*)'spettro di',levn(iv)%j2,(-1)**levn(iv)%l,'ev=', 
     &  levn(iv)%ei
       write(335,*)'spettro di',levn(iv)%j2,(-1)**levn(iv)%l,'ev=', 
     &  levn(iv)%ei
       do i=1,no+1
         write(334,*) i,wr(i),wi(i)
       end do
       allocate(wro(no+1))
       call order(wr,no+1,wro)
        do i=1,no+1
         write(335,*) i,wro(i)
       end do

       close(334)
       close(335)


       contains
        


      

!**********************************************************************
!**********************************************************************  
       subroutine permut(ndim,matr,mxt)

      implicit double precision (a-h,o-z)

      double precision, dimension(:,:), allocatable :: matr
      double precision, dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: mxt,nxt
      write(*,*) matr(1,1)
      !write(*,*) 'ndim',ndim,mxt(1)
      allocate (work(ndim))
      work=0.d0
      allocate(nxt(ndim))

      do i=1,ndim
       nxt(i)=i
      enddo

      do i=1,ndim
       work(:)=matr(i,:)
       do j=i,ndim
        if (nxt(j).eq.mxt(i)) then
        matr(i,:)=matr(j,:)
        matr(j,:)=work(:)
        nxt(j)=nxt(i)
        nxt(i)=mxt(i)
       endif
       enddo
      enddo

      do i=1,ndim
       nxt(i)=i
      enddo

      do i=1,ndim
       work(:)=matr(:,i)
       do j=i,ndim
        if (nxt(j).eq.mxt(i)) then
        matr(:,i)=matr(:,j)
        matr(:,j)=work(:)
        nxt(j)=nxt(i)
        nxt(i)=mxt(i)
        endif
       enddo
      enddo



      deallocate(work)

      return
      end subroutine
!******************************************************
!******************************************************
!******************************************************  
       subroutine permutv(ndim,matr,mxt)

      implicit double precision (a-h,o-z)

      double precision, dimension(:), allocatable :: matr
      double precision :: work

      integer, dimension(:), allocatable :: mxt,nxt
      write(*,*) matr(1)
      !write(*,*) 'ndim',ndim,mxt(1)
      !    allocate (work(ndim))
      work=0.d0
      allocate(nxt(ndim))

      do i=1,ndim
       nxt(i)=i
      enddo

      do i=1,ndim
       work=matr(i)
       do j=i,ndim
        if (nxt(j).eq.mxt(i)) then
        matr(i)=matr(j)
        matr(j)=work
        nxt(j)=nxt(i)
        nxt(i)=mxt(i)
       endif
       enddo
      enddo



!      deallocate(work)

      return
      end subroutine
!*********************************************************
      subroutine diagonalization(A,er,ei,v,i)

       include 'define.inc'

       double precision :: A(i,i),er(i),ei(i),v(i,i)
       double precision :: wor(max(1,6*i))

       icou1=max(1,6*i)
      er=0.d0
      ei=0.d0
      v=0.d0
      call DGEEV('N','V',i,A,max(1,i),er,ei,vl,
     &       max(1,i),v,max(1,i),wor,icou1,info)

       if(info.ne.0) then
        write(*,*) 'The diagonalization failed!!!'
        stop
       endif

       return
      end SUBROUTINE
!********************************************************************************
!********************************************************************************


        subroutine order(vec,mm,wro)

        integer ::mm,ipoz(mm) 
        double precision ::vec(mm)
        double precision, allocatable :: wro(:)
         !allocate(wro(mm))
        
        
      do i=1,mm

        ipoz(i)=i

        wro(i)=vec(i)

      enddo



      do i=1,mm

        do j=1,mm-i

          if (wro(j).gt.wro(j+1)) then

            xe=wro(j)

            wro(j)=wro(j+1)

            wro(j+1)=xe

            ipozz=ipoz(j)

            ipoz(j)=ipoz(j+1)

            ipoz(j+1)=ipozz

      endif

      enddo

      enddo
       return
      end subroutine  
!**************************************************************
!**************************************************************     
      end program
