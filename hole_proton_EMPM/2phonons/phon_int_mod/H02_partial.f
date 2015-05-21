
       Program H02_partial

!********************************************************************
!this program generate the interaction part of EMPM for odd nuclei
!********************************************************************

       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
 
       integer :: ii,jj,ipp1,inn1,imax1,ip
       Integer :: ipp2,inn2,imax2
       integer ::jalfa,lam,v,jv,Jlam,plam,ilam
       integer ::lam1,p,Jlam1,plam1,ilam1
       Integer :: p1,h1,jp1,jh1
       double precision :: s1,vv,Fint,racc
       double precision, allocatable :: rac(:,:,:,:,:)

       call input
       call Vinteraction

        jv=levp(iv)%j2

        do jh=0,jmax
         do jalfa=0,jjmax2
          do Jlam1=0,jjmax
           do Jlam=0,jjmax
            do jp=0,jmax
        rac(jh,jp,jalfa,Jlam,Jlam1)=racah(0.5d0*dble(jh),0.5d0*dble(jp),
     & dble(jalfa),dble(Jlam),dble(Jlam1),0.5d0*dble(jv))
            end do
           end do
          end do 
         end do
        end do



        open(22,file='int_par.dat',status='unknown',form='unformatted')
       Allocate(vinta(ihpmin:ihpmax,num_ph,num_ph,jjmin2:jjmax2))

         vinta=0.d0

       do jalfa=jjmin2,jjmax2

        do ih=ihpmin,ihpmax
         jp=levp(ip)%j2

       do ii=1,num_ph
       ilam=ph(ii)%lam
       Jlam=ph(ii)%Jph 
       plam=ph(ii)%par

       if(allocated(c1)) deallocate(c1)

       call cTDA(ilam,plam,Jlam)

       allocate(c1(1:ipnmax,1:ipnmax,-1:1)) 
       c1=0.d0
       c1=cph


       deallocate(cph)

       do ll=1,num_ph
       Jlam1=ph(ll)%Jph
       plam1=ph(ll)%par
       ilam1=ph(ll)%lam
        vv=0.d0
!**************************************************************************

         do pp=ippmin,ippmax
         jp=levn(pp)%j2
         racc=rac(jh,jp,jalfa,Jlam,Jlam1)
         if(dabs(rac).gt.0.d0) then
          vv=vv+c1(iv,hh,-1)*racc*Fcoup(pp,ih,ll)
          end if
        end do !hh

!**************************************************************************
         if(dabs(vv).gt.0.d0) then
         Vinta(ih,ii,ll,jalfa)=vv
         write(22) ih,ii,ll,jalfa,vv
         end if
        end do!ll
        end do!ii
        end do
        end do
    
        
        end program 
