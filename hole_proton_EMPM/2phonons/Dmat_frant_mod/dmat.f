       program Metricmat
 
       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc' 

       Integer :: ind1,ind2,j1,j2,sigma,il,ilj
       Integer :: ir,irj
       Integer :: lamp,lam,Jlp,Jl,ipp,ip,ipa,Jla,ilama
 
       double precision :: ep1,ep2
       double precision :: elp,el  
 
       double precision :: a1,a2,a3,a4,phase,dm
       double precision, allocatable :: rac(:,:,:,:,:)
       double precision,allocatable :: work(:),wr(:),wi(:),vr(:,:)
       integer :: info

       call input_files
       call coup_basis
   
       write(*,*) 'coup basis generated'
       Allocate(Dmat(dimbas,dimbas))

       Dmat=0.d0

       open(11,file='dmat.dat',status='unknown',form='formatted')
       open(33,file='dmateig.dat',status='unknown',form='formatted')


      
      Allocate(rac(0:jmax,0:jmax,0:jjmax2,0:jjmax2,0:2*jjmax2))
       rac=0.d0
      do i=0,jjmax!sigma
        do j=0,jmax
         do k=0,jjmax
          do l=0,jmax
           do m=0,jjmax
      rac(i,j,k,l,m)= racah(dble(i),0.5d0*dble(j),dble(k),0.5d0*dble(Jtot),
     & 0.5d0*dble(l),dble(m))
           end do
          end do
         end do
        end do
       end do


        do jj=1,dimbas
       write(*,*) jj,'of',dimbas 
       ilj=bas(jj)%ilam
       irj=bas(jj)%ir


       ind2=irj
       j2=levp(irj)%j2
       ep2=levp(irj)%ei
       lam=ph2(ilj)%alfa
       Jl=ph2(ilj)%J
       el=ph2(ilj)%el
       ip=ph2(ilj)%par

        call read_density(ilj,Jl,ip,irho)

       do ii=1,dimbas

       il=bas(ii)%ilam
       ir=bas(ii)%ir

       ind1=ir
       j1=levp(ir)%j2
       ep1=levp(ir)%ei
       lamp=ph2(il)%alfa
       Jlp=ph2(il)%J
       elp=ph2(il)%el
       ipp=ph2(il)%par
        if (ii.eq.jj) a1=1.d0
        if (ii.ne.jj) a1=0.d0
        dm=0.d0
       
      do sigma=abs(Jlp-Jl),(Jlp+Jl)

             a2=dsqrt(dble(2*sigma+1))
             a3=rac(sigma,j2,Jlp,j1,Jl)
             a4=den(il,ind2,ind1,sigma)
             dm=dm+a2*a3*phase*a4

      end do
      dm=a1+dm
      Dmat(ii,jj)=dm

        end do
       end do
      do il=1,dimbas
       do jl=1,dimbas
        write(11,*) il,jl,Dmat(il,jl)
          if(abs(Dmat(il,jl)-Dmat(jl,il)).ge.1.d-12)then
         write(3,*) il,jl,dabs(Dmat(il,jl)-Dmat(jl,il))
       end if    
      end do
      end do 

      allocate(vr(dimbas,dimbas))
      allocate (work(max(1,6*dimbas)))
      allocate(wr(dimbas),wi(dimbas))

      lwork=max(1,6*dimbas)
      call DSYEV('V','l',dimbas,Dmat,max(1,dimbas),wr,work,lwork,info)
       tr=0.d0
      do i=1,dimbas
       write(33,*) i,wr(i)!,wi(i)
      tr=tr+wr(i)
      end do
      write(*,*) 'trace of Dmatrix=',tr

      close(11)
      close(33)
!************************************************
      deallocate(levp,levn)
      deallocate(ph2)
      deallocate(wr,wi,vr)
!************************************************
      contains


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




      end program
