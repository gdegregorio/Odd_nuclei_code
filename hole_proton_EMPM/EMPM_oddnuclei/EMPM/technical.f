      MODULE technical

       include 'typedef.inc'

        PUBLIC :: levp, levn, lhfp, lhfn,
     &            Vpp,Vnn,Vpn,Fpp,Fnn,Fpn,
     &            cp,cn,eps,bas,dmat

        type(level_type),dimension(:),allocatable, save :: levp,levn
        type(ph_type),   dimension(:),allocatable, save :: ph
        type(dens_type) , dimension(:), allocatable, save :: den,h_den
        type(dimen_type),dimension(:),allocatable,save :: dime

        type(indbas_type), dimension(:),allocatable,save :: bas
         
        Integer, Save :: id,dimbas,no,iv

        Integer, save :: Jtot,ipartot 
        Integer, save :: partc,state  
        Integer, save :: pp1,pp2    

        Double precision,save ::eps

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpp(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fnn(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpn(:,:,:,:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Vcall(:,:)
       
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: dmat(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: d1(:,:)
 
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: AD(:,:)        
 
 
 
                
       END MODULE technical
