       Subroutine coup_basis
 
       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc' 
  
       integer :: jlam,jp,il,ip
       integer :: ii,jj
       integer :: dimb
     
        dimb=0
        dimb=id*dim2b+id
      write(*,*) 'dimb=',dimb      
      allocate(bas(dimb))
      write(*,*)'gen. of coup. basis for id=',id,'and num p=',dim2b
       dimbas=0
      open(1,file='infobase.dat',status='unknown',form='formatted')
       do ii=ippmin,ippmax
        do jj=1,num_2ph        
         jlam=ph2(jj)%J
         jp=levp(ii)%j2
         il=ph2(jj)%par
         ip=levp(ii)%ipar
         if(ip*il.eq.ipartot.and.ph2(jj)%el.lt.en_cut) then
          if(abs(2*jlam-jp).le.Jtot.and.Jtot.le.(2*jlam+jp)) then
           dimbas=dimbas+1
           bas(dimbas)%ilam=jj
           bas(dimbas)%ir=ii
          end if
         end if
        end do
       end do
       do i=1,dimbas
       write(1,*) i,bas(i)%ir,bas(i)%ilam
       end do
       close(1)

       write(*,*) 'basis of dimension =',dimbas,'generated'

       end subroutine
  
