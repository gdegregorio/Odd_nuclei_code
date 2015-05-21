      MODULE technical

       include 'typedef.inc'

        PUBLIC :: levp, levn, lhfp, lhfn,
     &            Vpp,Vnn,Vpn,Fpp,Fnn,Fpn,
     &            cp,cn,eps,bas,dmat,en_thr

        type(level_type),dimension(:),allocatable, save :: levp,levn
        !type(ph_type)   ,dimension(:),allocatable, save :: ph
        type(ph2_type)  ,dimension(:),allocatable, save :: ph2
        type(dens_type) ,dimension(:), allocatable, save :: den,h_den
        type(indbas_type), dimension(:),allocatable,save :: bas2ph
         type(dimen_type),dimension(:),allocatable,save ::dime

        Integer, Save :: id,dimbas,no,iv,dim2b!,dim1ph,dim2ph

        Integer, save :: Jtot,ipartot 
        Integer, save :: partc,state  
        Integer, save :: pp1,pp2    
   
        integer, allocatable,save ::good(:)
        Double precision,save ::eps,en_thr,Vcall

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpp(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fnn(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpn(:,:,:,:,:)

        DOUBLE PRECISION,allocatable,SAVE :: dime2ph(:,:)
       
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: dmat(:,:),DAD(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: d1(:,:),AD(:,:),d1m(:,:)
 
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: xx(:,:,:),C2ph(:,:)
 
 
 
                
       END MODULE technical
