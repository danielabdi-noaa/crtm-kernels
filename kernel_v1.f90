PROGRAM test_kernels

  IMPLICIT NONE

  INTEGER, PARAMETER :: fp = SELECTED_REAL_KIND(15)
  REAL(fp), PARAMETER :: ZERO          =  0.0_fp
  REAL(fp), PARAMETER :: ONE           =  1.0_fp
  REAL(fp), PARAMETER :: TWO           =  2.0_fp
  REAL(fp), PARAMETER :: OPTICAL_DEPTH_THRESHOLD = 0.000001_fp
  INTEGER, PARAMETER :: MAX_N_LAYERS   = 200
  INTEGER, PARAMETER :: INVALID_SENSOR = 0
  INTEGER, PARAMETER :: RT_ADA = 56
  INTEGER, PARAMETER :: MAX_N_ANGLES = 16
  INTEGER, PARAMETER :: MAX_N_LEGENDRE_TERMS = 16
  INTEGER, PARAMETER :: MAX_N_DOUBLING = 55
  INTEGER, PARAMETER :: MAX_N_SOI_ITERATIONS = 75
 
  !---- RTV type ---!
  TYPE :: RTV_type
  
    INTEGER :: n_Layers         = 0       ! Total number of atmospheric layers
    INTEGER :: n_Angles         = 0       ! Number of angles to be considered
    INTEGER :: n_SOI_Iterations = 0       ! Number of SOI iterations
    
    ! Planck radiances
    REAL(fp)                               :: Planck_Surface    = ZERO
    REAL(fp), DIMENSION(  0:MAX_N_LAYERS ) :: Planck_Atmosphere = ZERO

    ! Quadrature information
    REAL(fp), DIMENSION( MAX_N_ANGLES ) :: COS_Angle  = ZERO  ! Gaussian quadrature abscissa
    REAL(fp), DIMENSION( MAX_N_ANGLES ) :: COS_Weight = ZERO  ! Gaussian quadrature weights
    
    ! Scattering, visible model variables    
    INTEGER :: n_Streams         = 0       ! Number of *hemispheric* stream angles used in RT    

    !-----------------------------------
    ! Variables used in the ADA routines
    !-----------------------------------
    ! Flag to indicate the following arrays have all been allocated
    LOGICAL :: Is_Allocated = .FALSE.
     
    ! Phase function variables
    ! Forward and backward scattering phase matrices
    REAL(fp), ALLOCATABLE :: Pff(:,:,:)  ! MAX_N_ANGLES, MAX_N_ANGLES+1, MAX_N_LAYERS
    REAL(fp), ALLOCATABLE :: Pbb(:,:,:)  ! MAX_N_ANGLES, MAX_N_ANGLES+1, MAX_N_LAYERS

    !-----------------------------------
    ! Variables used in the SOI routines
    !-----------------------------------
    INTEGER :: Number_SOI_Iter = 0

    INTEGER , ALLOCATABLE :: Number_Doubling(:)  ! n_Layers
    REAL(fp), ALLOCATABLE :: Delta_Tau(:)        ! n_Layers
    REAL(fp), ALLOCATABLE :: Refl(:,:,:,:)       ! n_Angles, n_Angles, 0:MAX_N_DOUBLING, n_Layers
    REAL(fp), ALLOCATABLE :: Trans(:,:,:,:)      ! n_Angles, n_Angles, 0:MAX_N_DOUBLING, n_Layers
    REAL(fp), ALLOCATABLE :: Inv_BeT(:,:,:,:)    ! n_Angles, n_Angles, 0:MAX_N_DOUBLING, n_Layers
    REAL(fp), ALLOCATABLE :: C1(:,:)             ! n_Angles, n_Layers
    REAL(fp), ALLOCATABLE :: C2(:,:)             ! n_Angles, n_Layers

  END TYPE RTV_type

  !------- Test -------!
  TYPE(RTV_type) :: RTV
  INTEGER, PARAMETER :: N_LAYERS = MAX_N_LAYERS
  INTEGER :: k, alloc_stat
  INTEGER :: count_rate, count_start, count_end
  REAL :: elapsed

  !---- local arrays --- !
  REAL(fp), ALLOCATABLE, DIMENSION( :,:,: ) :: Pff_AD,Pbb_AD,s_Refl_AD,s_Trans_AD
  REAL(fp), ALLOCATABLE, DIMENSION( :,: ) :: s_source_UP_AD,s_source_DOWN_AD
  REAL(fp), DIMENSION(N_LAYERS) ::  w, T_OD
  REAL(fp), DIMENSION(N_LAYERS) :: w_AD, T_OD_AD
  REAL(fp), DIMENSION( 0:N_LAYERS ) :: Planck_Atmosphere_AD 
!$acc declare create(w, T_OD, w_AD, T_OD_AD, Planck_Atmosphere_AD)

  ALLOCATE(Pff_AD(MAX_N_ANGLES, MAX_N_ANGLES+1, N_LAYERS), &
           Pbb_AD(MAX_N_ANGLES, MAX_N_ANGLES+1, N_LAYERS), &
           s_Refl_AD(MAX_N_ANGLES, MAX_N_ANGLES, N_LAYERS), &
           s_Trans_AD(MAX_N_ANGLES, MAX_N_ANGLES, N_LAYERS), &
           s_source_UP_AD(MAX_N_ANGLES, N_LAYERS), &
           s_source_DOWN_AD(MAX_N_ANGLES, N_LAYERS), &
           STAT = alloc_stat )
  IF ( alloc_stat /= 0 ) STOP

  !---- fill arrays with random numbers ----!
  PRINT*, "Filling arrays with random values"
  CALL RANDOM_NUMBER(Pff_AD)
  CALL RANDOM_NUMBER(Pbb_AD)
  CALL RANDOM_NUMBER(s_Refl_AD)
  CALL RANDOM_NUMBER(s_Trans_AD)
  CALL RANDOM_NUMBER(s_source_UP_AD)
  CALL RANDOM_NUMBER(s_source_DOWN_AD)
  CALL RANDOM_NUMBER(w)
  CALL RANDOM_NUMBER(T_OD)
  CALL RANDOM_NUMBER(w_AD)
  CALL RANDOM_NUMBER(T_OD_AD)
  CALL RANDOM_NUMBER(Planck_Atmosphere_AD)
  PRINT*, "Finished filling arrays with random values."

!$acc update device(w, T_OD, w_AD, T_OD_AD, Planck_Atmosphere_AD)
!$acc enter data copyin(Pff_AD,Pbb_AD,s_Refl_AD,s_Trans_AD, &
!$acc                  s_source_UP_AD,s_source_DOWN_AD)

  !---- create RTV array ----!
  PRINT*, "Creating RTV"
  RTV%n_Streams = 6
  CALL RTV_Create( RTV, MAX_N_ANGLES, MAX_N_LEGENDRE_TERMS, N_LAYERS )
  PRINT*, "Finished creating RTV"


  !---- call kernel ----!
  PRINT*, "Calling kernel"  
  CALL SYSTEM_CLOCK (count_rate=count_rate)
  CALL SYSTEM_CLOCK (count=count_start)

  DO k = 1, N_LAYERS
  CALL CRTM_Doubling_layer_AD(RTV%n_Streams, RTV%n_Angles, k, w( k ), T_OD( k ),      &        !Input
                           RTV%COS_Angle, RTV%COS_Weight, RTV%Pff( :, :, k ), RTV%Pbb( :, :, k ), & ! Input
                           RTV%Planck_Atmosphere( k ),    & !Input
                           s_trans_AD( :, :, k ), s_refl_AD( :, :, k ), s_source_up_AD( :, k ),   & 
                           s_source_down_AD( :, k ), RTV, w_AD( k ), T_OD_AD( k ), Pff_AD( :, :, k ), & 
                           Pbb_AD( :, :, k ), Planck_Atmosphere_AD( k ) )  !Output
  ENDDO

  CALL SYSTEM_CLOCK (count=count_end)
  elapsed = REAL (count_end - count_start) / REAL (count_rate)
  PRINT*
  PRINT*
  PRINT*, "Finished executing kernel in =", elapsed  
  PRINT*
  !------- Test -------!
!$acc update self(s_trans_AD)
  PRINT*, s_trans_AD(:,1,1)

 CONTAINS

  SUBROUTINE  myMATMUL(A, B, C)
    REAL(fp), INTENT(IN), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES) :: A, B
    REAL(fp), INTENT(OUT), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES) :: C
!$acc declare copyin(A,B) copyout(C)
    REAL(fp) :: acc
    INTEGER :: I, J, K

!$acc parallel
!$acc loop collapse(2) private(acc)
    DO I = 1, MAX_N_ANGLES
      DO J = 1, MAX_N_ANGLES
        acc = 0
!$acc loop reduction(+:acc)
        DO K = 1, MAX_N_ANGLES
            acc = acc + A(I,K) * B(K,J)
        END DO
        C(I, J) = acc
      END DO
    END DO
!$acc end parallel

  END SUBROUTINE myMATMUL

  SUBROUTINE myTRANSPOSE(A, At)
    REAL(fp), INTENT(IN), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES) :: A
    REAL(fp), INTENT(OUT), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES) :: At
    INTEGER :: I, J
!$acc declare copyin(A) copyout(At)

!$acc parallel
!$acc loop collapse(2)
    DO I = 1, MAX_N_ANGLES
      DO J = 1, MAX_N_ANGLES
         At(I,J) = A(J,I)
      END DO
    END DO
!$acc end parallel

  END SUBROUTINE myTRANSPOSE

  SUBROUTINE RTV_Create( &
    RTV, &
    n_Angles        , &
    n_Legendre_Terms, &
    n_Layers          )
    ! Arguments
    TYPE(RTV_type), INTENT(OUT) :: RTV
    INTEGER       , INTENT(IN)  :: n_Angles        
    INTEGER       , INTENT(IN)  :: n_Legendre_Terms
    INTEGER       , INTENT(IN)  :: n_Layers        
    ! Local variables
    INTEGER :: alloc_stat

    ! Check input
    IF ( n_Angles < 1 .OR. n_Legendre_Terms < 1 .OR. n_Layers < 1 ) RETURN
    
    ALLOCATE( RTV%Pff(n_Angles, n_Angles+1, n_Layers) , &
         RTV%Pbb(n_Angles, n_Angles+1, n_Layers) , &
         STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN

    ! Perform the allocation for SOI variables
    ALLOCATE( RTV%Number_Doubling(n_Layers), &
              RTV%Delta_Tau(n_Layers), &      
              RTV%Refl(n_Angles, n_Angles, 0:MAX_N_DOUBLING, n_Layers), &
              RTV%Trans(n_Angles, n_Angles, 0:MAX_N_DOUBLING, n_Layers), &
              RTV%Inv_BeT(n_Angles, n_Angles, 0:MAX_N_DOUBLING, n_Layers), &
              RTV%C1(n_Angles, n_Layers), &
              RTV%C2(n_Angles, n_Layers), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN

    ! Set dimensions
    RTV%n_Layers         = n_Layers
    RTV%n_Angles         = n_Angles
    
    RTV%n_SOI_Iterations = 0

    ! Set the allocate flag
    RTV%Is_Allocated = .TRUE.

    !------ Fill arrays with random values --- !
    CALL RANDOM_NUMBER(RTV%Pff)
    CALL RANDOM_NUMBER(RTV%Pbb)
    RTV%Number_Doubling(:) = MAX_N_DOUBLING
    CALL RANDOM_NUMBER(RTV%Delta_Tau)
    CALL RANDOM_NUMBER(RTV%Refl)
    CALL RANDOM_NUMBER(RTV%Trans)
    CALL RANDOM_NUMBER(RTV%Inv_BeT)
    CALL RANDOM_NUMBER(RTV%COS_Angle)
    CALL RANDOM_NUMBER(RTV%COS_Angle)
    CALL RANDOM_NUMBER(RTV%C1)
    CALL RANDOM_NUMBER(RTV%C2)
    RTV%Refl=RTV%Refl/100.0_fp
    RTV%Trans=RTV%Trans/100.0_fp
    RTV%Inv_Bet=RTV%Inv_Bet/100.0_fp
    !------  Deep copy -------! 
!$acc enter data copyin(RTV)
!$acc enter data copyin(RTV%Pff,RTV%Pbb,RTV%Number_Doubling,RTV%Delta_Tau, &
!$acc                RTV%Refl,RTV%Trans,RTV%Inv_BeT,RTV%C1,RTV%C2)

  END SUBROUTINE RTV_Create

  SUBROUTINE CRTM_Doubling_layer_AD(n_streams, & ! Input, number of streams
                                         NANG, & ! Input, number of angles
                                           KL, & ! Input, number of angles
                                single_albedo, & ! Input, single scattering albedo
                                optical_depth, & ! Input, layer optical depth
                                    COS_Angle, & ! Input, COSINE of ANGLES
                                   COS_Weight, & ! Input, GAUSSIAN Weights
                                           ff, & ! Input, Phase matrix (forward part)
                                           bb, & ! Input, Phase matrix (backward part)
                                  Planck_Func, & ! Input, Planck for layer temperature
                                     trans_AD, & ! Input, layer tangent-linear trans 
                                      refl_AD, & ! Input, layer tangent-linear refl 
                                 source_up_AD, & ! Input, layer tangent-linear source_up 
                               source_down_AD, & ! Input, layer tangent-linear source_down 
                                          RTV, & ! Input, structure containing forward results 
                             single_albedo_AD, & ! Output adjoint single scattering albedo
                             optical_depth_AD, & ! Output AD layer optical depth
                                        ff_AD, & ! Output AD forward Phase matrix
                                        bb_AD, & ! Output AD backward Phase matrix
                               Planck_Func_AD)   ! Output AD Planck for layer temperature


    INTEGER, INTENT(IN) :: n_streams,NANG,KL
    TYPE(RTV_type), INTENT(IN) :: RTV
!$acc declare copyin(RTV)
    REAL(fp), INTENT(IN), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES+1) :: ff,bb
    REAL(fp), INTENT(IN), DIMENSION(1:MAX_N_ANGLES) :: COS_Angle, COS_Weight 
    REAL(fp), INTENT(IN) :: single_albedo,optical_depth,Planck_Func
!$acc declare copyin(ff,bb, COS_Angle, COS_Weight)

    ! Tangent-Linear Part
    REAL(fp), INTENT( INOUT ), DIMENSION( 1:MAX_N_ANGLES,1:MAX_N_ANGLES ) :: trans_AD,refl_AD
    REAL(fp), INTENT( INOUT ), DIMENSION( 1:MAX_N_ANGLES ) :: source_up_AD,source_down_AD
    REAL(fp), INTENT( INOUT ) :: single_albedo_AD
    REAL(fp), INTENT( INOUT ) :: optical_depth_AD,Planck_Func_AD
    REAL(fp), INTENT(INOUT), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES+1) :: ff_AD,bb_AD
!$acc declare copy(trans_AD,refl_AD,source_up_AD,source_down_AD,ff_AD,bb_AD)

    ! internal variables
    REAL(fp), DIMENSION(MAX_N_ANGLES,MAX_N_ANGLES) :: term1,term2,term3,term4,term5_AD
    REAL(fp), DIMENSION(MAX_N_ANGLES,MAX_N_ANGLES) :: trans1,trans3,trans4,temp1,temp2,temp3
!$acc declare create(term1,term2,term3,term4,term5_AD,trans1,trans3,trans4,temp1,temp2,temp3)
    REAL(fp) :: s, c
    REAL(fp) :: s_AD, c_AD, Delta_Tau_AD
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: C1_AD, C2_AD
    INTEGER :: i,j,L
!$acc declare create(C1_AD,C2_AD)

    ! Tangent-Linear Beginning
    IF( optical_depth < OPTICAL_DEPTH_THRESHOLD ) THEN
!$acc kernels
      trans_AD = ZERO
      refl_AD = ZERO
      source_up_AD = ZERO
      source_down_AD = ZERO
!$acc end kernels
      RETURN
    ENDIF
   
!$acc kernels
!$acc loop reduction(+:Planck_Func_AD)
    DO i = NANG, 1, -1
      source_up_AD(i) = source_up_AD(i) + source_down_AD(i)
      source_down_AD(i) = ZERO
      C2_AD(i) = -source_up_AD(i)*Planck_Func
      C1_AD(i) = -source_up_AD(i)*Planck_Func
      Planck_Func_AD = Planck_Func_AD + (ONE-RTV%C1(i,KL)-RTV%C2(i,KL))*source_up_AD(i)
    END DO

    ! Compute the source function in the up and downward directions.
    IF(NANG == (n_Streams+1)) THEN
        trans_AD(NANG,NANG)=trans_AD(NANG,NANG)+C1_AD(NANG)
    ENDIF
!$acc loop collapse(2)
    DO i = NANG, 1, -1
      DO j = n_Streams, 1, -1 
        refl_AD(i,j)=refl_AD(i,j)+C2_AD(i)
        trans_AD(i,j)=trans_AD(i,j)+C1_AD(i)
      END DO
    END DO
!$acc end kernels

    DO L = RTV%Number_Doubling(KL), 1, -1 

      CALL myMATMUL(RTV%Trans(:,:,L-1,KL),RTV%Inv_BeT(:,:,L,KL), term1)
      CALL myMATMUL(RTV%Inv_BeT(:,:,L,KL),RTV%Refl(:,:,L-1,KL), term2)
      CALL myMATMUL(RTV%Inv_BeT(:,:,L,KL),RTV%Trans(:,:,L-1,KL), term3)
      CALL myMATMUL(term2,RTV%Trans(:,:,L-1,KL), term4)

      CALL myTRANSPOSE(term1,trans1)
      CALL myTRANSPOSE(term3,trans3)
      CALL myTRANSPOSE(term4,trans4)

      CALL myMATMUL(trans1,trans_AD,temp1)
      CALL myMATMUL(temp1,trans3,term5_AD)

      CALL myMATMUL(trans_AD,trans3,temp1)
      CALL myMATMUL(trans1,trans_AD,temp2)
!$acc kernels
      trans_AD = temp1 + temp2
!$acc end kernels

      CALL myMATMUL(term1,RTV%Refl(:,:,L-1,KL),temp1)
      CALL myTRANSPOSE(temp1,temp2)
      CALL myMATMUL(temp2,refl_AD,temp1) 
!$acc kernels
      trans_AD = trans_AD + temp1
!$acc end kernels

      CALL myMATMUL(trans1,refl_AD,temp2)
      CALL myMATMUL(temp2,trans4,temp1)
!$acc kernels
      term5_AD = term5_AD + temp1
!$acc end kernels

      CALL myMATMUL(refl_AD,trans4,temp1)
!$acc kernels
      trans_AD = trans_AD + temp1
!$acc end kernels

      CALL myMATMUL(trans1,refl_AD,temp1)
      CALL myTRANSPOSE(RTV%Trans(:,:,L-1,KL),temp3)
      CALL myMATMUL(temp1,temp3,temp2)
!$acc kernels
      refl_AD = refl_AD + temp2
!$acc end kernels

      CALL myTRANSPOSE(RTV%Refl(:,:,L-1,KL),temp3)
      CALL myMATMUL(term5_AD,temp3,temp1)
!$acc kernels
      refl_AD = refl_AD + temp1
!$acc end kernels

      CALL myMATMUL(temp3,term5_AD, temp1)
!$acc kernels
      refl_AD = refl_AD + temp1
!$acc end kernels

    ENDDO

    s = RTV%Delta_Tau(KL) * single_albedo
    c_AD = ZERO
    s_AD = ZERO
    Delta_Tau_AD=ZERO

!$acc parallel loop firstprivate(c_AD) reduction(+:Delta_Tau_AD,s_AD)
    DO i = NANG, 1, -1

      c = s/COS_Angle(i)
      Delta_Tau_AD = Delta_Tau_AD - trans_AD(i,i)/COS_Angle(i)

!$acc loop independent reduction(+:c_AD)
      DO j = NANG, 1, -1
        c_AD = c_AD + trans_AD(i,j)*ff(i,j)*COS_Weight(j)
        ff_AD(i,j)=ff_AD(i,j)+trans_AD(i,j)*c*COS_Weight(j)
        c_AD = c_AD + refl_AD(i,j)*bb(i,j)*COS_Weight(j)
        bb_AD(i,j)=bb_AD(i,j) + refl_AD(i,j)*c*COS_Weight(j)
      END DO

      s_AD = s_AD + c_AD/COS_Angle(i) 
      c_AD = ZERO

    ENDDO
!$acc end parallel

    Delta_Tau_AD = Delta_Tau_AD + s_AD* single_albedo
    single_albedo_AD = single_albedo_AD+RTV%Delta_Tau(KL) * s_AD
    optical_depth_AD = optical_depth_AD + Delta_Tau_AD/(TWO**RTV%Number_Doubling(KL))

  END SUBROUTINE CRTM_Doubling_layer_AD

END PROGRAM test_kernels