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
      call xmatr
      call permutv(dimbas,Fcoup,mxt)

      allocate(term(no))

       term=0.d0
       open(897,file='coup_term.dat',status='unknown',form='formatted')

       do j=1,no
        vec=0.d0
 
        do io=1,dimbas
         i=mxt1(io) 
         Jl=ph(bas(i)%ilam)%Jph
         jp=levp(bas(i)%ir)%j2
         phase=(-1)**(0.5d0*(Jtot)-0.5d0*jp-Jl)
         vec=vec+(dble(2*Jl+1))**0.5d0*phase*Fcoup(i)*xm(i,j)
        end do
        term(j)=(1.d0/(dble(levp(iv)%j2)+1))*vec
        write(897,*) j,term(j)
      
        end do
         close(897)

        allocate(ham(no+1,no+1))
         ham=0.d0
        ham(1,1)=levp(iv)%ei
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
       write(334,*)'spettro di',levp(iv)%j2,(-1)**levp(iv)%l,'ev=', 
     &  levp(iv)%ei
       write(335,*)'spettro di',levp(iv)%j2,(-1)**levp(iv)%l,'ev=', 
     &  levp(iv)%ei
       do i=1,no+1
         write(334,*) i,wr(i),wi(i)
        do j=1,no+1
        write(33,*) i,j,vr(j,i)
       end do
        end do
      open(22,file='struct_coup.dat',status='unknown',form='formatted')

      open(222,file='toteigv.dat',status='unknown',form='formatted')
       write(222,*) no+1
       do ii=1,no+1
       if(dabs(wr(ii)).lt.30.d0) then
        write(22,*)'***************************'
        write(22,*) ii,wr(ii)
        write(22,*)'***************************'
        write(22,*)'v,c,c^2'
        write(22,*) 1,vr(1,ii),vr(1,ii)**2
        !write(222,*) 1,ii,vr(1,ii)
         do j=2,no+1
        !write(222,*) j,ii,vr(j,ii)       
          if(dabs(vr(j,ii)).gt.0.01)then
           write(22,*) j,vr(j,ii),vr(j,ii)**2
          end if
         end do    
        end if
        end do
       do ii=1,no+1
        write(222,*) 1,ii,vr(1,ii)
         do j=2,no+1
          write(222,*) j,ii,vr(j,ii)  
         end do
       end do
    
       allocate(wro(no+1))
       call order(wr,no+1,wro)
        do i=1,no+1
         write(335,*) i,wro(i)
       end do
       close(222)
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
