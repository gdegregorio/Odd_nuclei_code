      MODULE technical

       include 'typedef.inc'

        PUBLIC :: levp, levn, lhfp, lhfn,
     &            Vpp,Vnn,Vpn,Fpp,Fnn,Fpn,
     &            cp,cn,eps,bas,dmat

        type(level_type),dimension(:),allocatable, save :: levp,levn
        type(ph_type),   dimension(:),allocatable, save :: ph
        type(dens_type) , dimension(:), allocatable, save :: den,h_den
        type(dimen_type),dimension(:),allocatable,save :: dime
        type(cph_type),dimension(:),allocatable,save :: cp,cn,C1p,C1n
        type(indbas_type), dimension(:),allocatable,save :: bas
         
        Integer, Save :: id,dimbas,no,iv

        Integer, save :: Jtot,ipartot 
        Integer, save :: partc,state  
        Integer, save :: pp1,pp2    
        Integer, allocatable,save ::mxt(:),mxt1(:)
        Double precision,save ::eps

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpp(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vnn(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vpn(:,:,:,:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpp(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fnn(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpn(:,:,:,:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fcoup(:),xm(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: dmat(:,:),d1(:,:)

         DOUBLE PRECISION,ALLOCATABLE,SAVE :: eigen(:)      
 
 
 
                
       END MODULE technical
