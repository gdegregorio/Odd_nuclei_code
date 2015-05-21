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

       Allocate(rac(0:2*jjmax,0:jmax,0:jjmax,0:jmax,0:jjmax))
       rac=0.d0
      do i=0,jjmax!sigma
        do j=0,jmax!h
         do k=0,jjmax!lam'
          do l=0,jmax!h'
           do m=0,jjmax!lam
      rac(i,j,k,l,m)= racah(dble(i),0.5d0*dble(j),dble(k),
     & 0.5d0*dble(Jtot),0.5d0*dble(l),dble(m))
        ! if(i.eq.1.and.j.eq.3.and.l.eq.3) then
        !  write(21,*) k,m,rac(i,j,k,l,m)
        ! end if
          end do
          end do
         end do
        end do
       end do

       do ii=1,dimbas

       write(*,*) ii,'of',dimbas

       ind1=indbase(ii)%ir!h'
       j1=levp(ind1)%j2
       lampp=indbase(ii)%ilam!lamp
       Jlp=ph(lampp)%Jph
       elp=ph(lampp)%el
       ipp=ph(lampp)%par
       lamp=ph(lampp)%lam

       call read_density(-1,-1,ipp,Jlp,lamp)

        do jj=1,dimbas
       ind2=indbase(jj)%ir!h
       j2=levp(ind2)%j2
       lamm=indbase(jj)%ilam!lambda
       Jl=ph(lamm)%Jph
       ip=ph(lamm)%par
       lam=ph(lamm)%lam
 
        if (ii.eq.jj) a1=1.d0
        if (ii.ne.jj) a1=0.d0
        dm=0.d0
          do sigma=abs(Jlp-Jl),(Jlp+Jl)
             a2=dsqrt(dble(2*sigma+1))
             a3=rac(sigma,j2,Jlp,j1,Jl)
             a4=den(lamm,ind1,ind2,sigma)
             dm=dm+a2*a3*a4
           ! if(ii.eq.1.and.jj.eq.46) then
           !  write(22,*) ii,jj,a1,a2,a3,a4,dm
           !  write(23,*) sigma,j2,Jlp,j1,Jl,a3
            !end if
           ! else if(ii.eq.46.and.jj.eq.1) then
           !write(22,*) ii,jj,a1,a2,a3,a4,dm
           ! write(23,*) sigma,j2,Jlp,j1,Jl,a3
           ! end if
          end do
      if(ii.eq.1.and.jj.eq.46) then 
      write(24,*) a1,dm
       end if
      dm=a1+dm
      !if(ii.eq.1.and.jj.eq.46) then
      
      !Dmat(ii,jj)=dm
      !if(ii.eq.1) then
      !write(22,*) Dmat(ii,jj)
      !end if  
      !Dmat(jj,ii)=dm
      write(11,*) ii,jj,dm!Dmat(ii,jj)
        end do
       end do
       !write(*,*) Dmat(1,46),Dmat(46,1)
       close(11)
      open(11,file='dmat.dat',status='unknown',form='formatted')
      do while(.not.eof(11)) 
       read(11,*) il,jl,Dmat(il,jl)
      end do
      close(11)
      open(40,file='sym_dmat.dat',status='unknown',form='formatted')
       do i=1,dimbas
        do j=1,dimbas
         if(dabs(Dmat(i,j)-Dmat(j,i)).gt.0.000000001) then
          write(40,*) 'nsym',i,j,Dmat(i,j),Dmat(j,i)
         end if
        end do 
       end do
      close(40)
      open(43,file='dmat_n0.dat',status='unknown',form='formatted')
      do il=1,dimbas
       do jl=1,dimbas
        !write(11,*) il,jl,Dmat(il,jl)
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

      !close(11)
      close(33)
!************************************************
      deallocate(indbase)
      deallocate(work,w,ph)
!************************************************

      end program
