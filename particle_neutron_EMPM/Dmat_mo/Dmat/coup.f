       Subroutine coup_basis
 
       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc' 
  
       integer :: jlam,jp,il,ip
       integer :: ii,jj,kk
       double precision ::en_th

      write(*,*)'gen. of coup. basis for id=',id,'and num p=',num_ph
      open(1,file='infobase.dat',status='unknown',form='formatted')
      write(*,*) 'energy threshold for 1 phonon states?'
      read(*,*) en_th
      write(*,*) 'energy='
       ind=0
       do kk=1,num_ph
        if(ph(kk)%el.le.en_th) then
         ind=ind+1
        end if
       end do

       do ii=ipnmin,ipnmax
        do jj=1,num_ph        
         jlam=ph(jj)%Jph
         jp=levn(ii)%j2
         il=ph(jj)%par
         ip=levn(ii)%ipar
         if(ip*il.eq.ipartot.and.ph(jj)%el.le.en_th) then
          if(abs(2*jlam-jp).le.Jtot.and.Jtot.le.(2*jlam+jp)) then
           dimbas=dimbas+1
          endif
         end if
        end do
       end do

       allocate(indbase(dimbas))
       kk=1
       do ii=ipnmin,ipnmax
       do jj=1,num_ph
       jlam=ph(jj)%Jph
       jp=levn(ii)%j2
       il=ph(jj)%par
       ip=levn(ii)%ipar
       if(ip*il.eq.ipartot.and.ph(jj)%el.le.en_th) then
       if(abs(2*jlam-jp).le.Jtot.and.Jtot.le.(2*jlam+jp)) then
         indbase(kk)%ir=ii
         indbase(kk)%ilam=jj
       write(1,*)kk,ii,jj
       kk=kk+1
          end if
         end if
        end do
       end do
       close(1)
       write(*,*) 'basis of dimension =',dimbas,'generated'
       !deallocate(levp,levn)
       deallocate(dime)
         !stop
        return
       end subroutine
  
