      MODULE technical

       include 'typedef.inc'

        PUBLIC :: levp, levn, lhfp, lhfn,
     &            Fpp,Fnn,Fpn,dime,dime2ph,
     &            cp,cn,eps,bas,dmat,en_thr

        type(level_type),dimension(:),allocatable, save :: levp,levn
        type(ph_type)   ,dimension(:),allocatable, save :: ph
        type(dens_type) ,dimension(:), allocatable, save :: den,h_den
        type(indbas_type), dimension(:),allocatable,save :: bas
        type(dimen_type),dimension(:),allocatable,save ::dime


        Integer, Save :: id,dimbas,no,iv,dim1ph,dimbs
         
        Integer, save :: Jtot,ipartot 
        Integer, save :: partc,state  
        Integer, save :: pp1,pp2    

        Double precision,save ::eps,en_thr,Vcall

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpp(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fnn(:,:,:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fpn(:,:,:,:,:)

        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Fcoup(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: vinta(:,:,:,:)

        double precision,allocatable, save ::eigv(:,:) 
        double precision,allocatable, save ::cph(:,:,:),c1(:,:,:)
        double precision,allocatable, save :: c2(:,:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: xx(:,:)
 
 
 
                
       END MODULE technical
