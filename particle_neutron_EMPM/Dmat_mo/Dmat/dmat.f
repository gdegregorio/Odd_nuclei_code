       program Metricmat
 
       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc' 

       Integer :: ind1,ind2,j1,j2,sigma,lampp,lamm
       Integer :: lamp,lam,Jlp,Jl,ipp,ip,ipa,Jla,ilama
 
       double precision :: ep1,ep2
       double precision :: elp,el  
 
       double precision :: a1,a2,a3,a4,phase,dm
       double precision,allocatable :: rac(:,:,:,:,:)
       double precision,allocatable :: work(:),w(:)
       integer :: info

       call input_files
       call coup_basis
   
       write(*,*) 'coup basis generated'
       Allocate(Dmat(dimbas,dimbas))

       Dmat=0.d0

       open(11,file='dmat.dat',status='unknown',form='formatted')
       open(33,file='dmateig.dat',status='unknown',form='formatted')

       Allocate(rac(0:jmax,0:jmax,0:jjmax,0:jjmax,0:2*jjmax))
       rac=0.d0
      do i=0,jmax
        do j=0,jmax
         do k=0,jjmax
          do l=0,jjmax
           do m=0,jjmax
      rac(i,j,k,l,m)= racah(0.5d0*dble(i),0.5d0*dble(j),dble(k),dble(l),
     &       dble(m),0.5d0*dble(Jtot))
           end do
          end do
         end do
        end do
       end do
  



       do ii=1,dimbas

       write(*,*) ii,'of',dimbas

       ind1=indbase(ii)%ir
       j1=levn(ind1)%j2
       lampp=indbase(ii)%ilam
       Jlp=ph(lampp)%Jph
       elp=ph(lampp)%el
       ipp=ph(lampp)%par
       lamp=ph(lampp)%lam
       call read_density(1,1,ipp,Jlp,lamp)
         !if(ii.eq.1224) then
         ! aojifa=0.d0
         !end if      
        do jj=ii,dimbas
         !if(ii.eq.1224) write(*,*) jj 
       ind2=indbase(jj)%ir
       j2=levn(ind2)%j2
       lamm=indbase(jj)%ilam
       Jl=ph(lamm)%Jph
       ip=ph(lamm)%par
       lam=ph(lamm)%lam
 
        if (ii.eq.jj) a1=1.d0
        if (ii.ne.jj) a1=0.d0
        dm=0.d0
       
          do sigma=abs(Jlp-Jl),(Jlp+Jl)

             a2=dsqrt(dble(2*sigma+1))
             a3=rac(j1,j2,Jlp,Jl,sigma)
!ah(0.5d0*dble(j1),0.5d0*dble(j2),dble(Jlp),dble(Jl),
!a&                dble(sigma),0.5d0*dble(Jtot))
             a4=den(lamm,ind2,ind1,sigma)
             phase=dfloat((-1)**((j2-Jtot)/2+Jl))
             dm=dm+a2*a3*phase*a4
          end do
      dm=a1-dm        
      Dmat(ii,jj)=dm
      Dmat(jj,ii)=dm
        end do
       end do
      open(43,file='dmat_n0.dat',status='unknown',form='formatted')
      do il=1,dimbas
       do jl=1,dimbas
        write(11,*) il,jl,Dmat(il,jl)
        if(Dmat(il,jl).ge.0.2) then
        write(43,*) il,jl,Dmat(il,jl)
        end if 
      end do
      end do 
      close(43)
      allocate(work(max(1,6*dimbas)))
      allocate(w(dimbas))

      lwork=max(1,6*dimbas)
      call DSYEV('V','l',dimbas,Dmat,max(1,dimbas),w,work,lwork,info)
       tr=0.d0
      do i=1,dimbas
       write(33,*) i,w(i)
      tr=tr+w(i)
      end do
      write(*,*) 'trace of Dmatrix=',tr

      close(11)
      close(33)
!************************************************
      deallocate(indbase)
      deallocate(work,w,ph)
!************************************************

      end program
