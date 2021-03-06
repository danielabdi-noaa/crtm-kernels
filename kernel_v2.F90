!----------------
! kernels
!----------------
MODULE my_kernels

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

  INTEGER :: N_GPUS 

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

 CONTAINS

  SUBROUTINE  myMATMUL(A, B, C)
!$acc routine vector
    REAL(fp), INTENT(IN), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES) :: A, B
    REAL(fp), INTENT(OUT), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES) :: C
    REAL(fp) :: acc
    INTEGER :: I, J, K

!$acc loop collapse(2) private(acc)
    DO I = 1, MAX_N_ANGLES
      DO J = 1, MAX_N_ANGLES
        acc = 0
!$acc loop seq
        DO K = 1, MAX_N_ANGLES
            acc = acc + A(I,K) * B(K,J)
        END DO
        C(I, J) = acc
      END DO
    END DO

  END SUBROUTINE myMATMUL

  SUBROUTINE myTRANSPOSE(A, At)
!$acc routine vector
    REAL(fp), INTENT(IN), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES) :: A
    REAL(fp), INTENT(OUT), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES) :: At
    INTEGER :: I, J

!$acc loop collapse(2)
    DO I = 1, MAX_N_ANGLES
      DO J = 1, MAX_N_ANGLES
         At(I,J) = A(J,I)
      END DO
    END DO

  END SUBROUTINE myTRANSPOSE

  SUBROUTINE myADD(A, B, C)
!$acc routine vector
    REAL(fp), INTENT(IN), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES) :: A, B
    REAL(fp), INTENT(OUT), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES) :: C
    INTEGER :: I, J

!$acc loop collapse(2)
    DO I = 1, MAX_N_ANGLES
      DO J = 1, MAX_N_ANGLES
         C(I,J) = A(I,J) + B(I,J)
      END DO
    END DO

  END SUBROUTINE myADD

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
    CALL RANDOM_NUMBER(RTV%COS_Weight)
    CALL RANDOM_NUMBER(RTV%C1)
    CALL RANDOM_NUMBER(RTV%C2)
    RTV%Refl=RTV%Refl/100.0_fp
    RTV%Trans=RTV%Trans/100.0_fp
    RTV%Inv_Bet=RTV%Inv_Bet/100.0_fp

    !------  Manual deep copy -------! 
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
                               Planck_Func_AD, & ! Output AD Planck for layer temperature
                               term1,term2,term3,term4,term5_AD,trans1,trans3,trans4,temp1,temp2,temp3,C1_AD,C2_AD)   ! Temporaries

!$acc routine vector
    INTEGER, INTENT(IN) :: n_streams,NANG,KL
    TYPE(RTV_type), INTENT(IN) :: RTV
    REAL(fp), INTENT(IN), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES+1) :: ff,bb
    REAL(fp), INTENT(IN), DIMENSION(1:MAX_N_ANGLES) :: COS_Angle, COS_Weight 
    REAL(fp), INTENT(IN) :: single_albedo,optical_depth,Planck_Func

    ! Tangent-Linear Part
    REAL(fp), INTENT( INOUT ), DIMENSION( 1:MAX_N_ANGLES,1:MAX_N_ANGLES ) :: trans_AD,refl_AD
    REAL(fp), INTENT( INOUT ), DIMENSION( 1:MAX_N_ANGLES ) :: source_up_AD,source_down_AD
    REAL(fp), INTENT( INOUT ) :: single_albedo_AD
    REAL(fp), INTENT( INOUT ) :: optical_depth_AD,Planck_Func_AD
    REAL(fp), INTENT(INOUT), DIMENSION(1:MAX_N_ANGLES,1:MAX_N_ANGLES+1) :: ff_AD,bb_AD

    ! internal variables
    REAL(fp), INTENT(OUT),DIMENSION(MAX_N_ANGLES,MAX_N_ANGLES) :: term1,term2,term3,term4,term5_AD
    REAL(fp), INTENT(OUT),DIMENSION(MAX_N_ANGLES,MAX_N_ANGLES) :: trans1,trans3,trans4,temp1,temp2,temp3
    REAL(fp), INTENT(OUT),DIMENSION(MAX_N_ANGLES) :: C1_AD, C2_AD
    REAL(fp) :: s, c
    REAL(fp) :: s_AD, c_AD, Delta_Tau_AD
    INTEGER :: i,j,L

    ! Tangent-Linear Beginning
    IF( optical_depth < OPTICAL_DEPTH_THRESHOLD ) THEN
      trans_AD = ZERO
      refl_AD = ZERO
      source_up_AD = ZERO
      source_down_AD = ZERO
      RETURN
    ENDIF
   
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

!$acc loop seq
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
      CALL myADD(temp1,temp2,trans_AD)

      CALL myMATMUL(term1,RTV%Refl(:,:,L-1,KL),temp1)
      CALL myTRANSPOSE(temp1,temp2)
      CALL myMATMUL(temp2,refl_AD,temp1) 
      CALL myADD(trans_AD,temp1,trans_AD)

      CALL myMATMUL(trans1,refl_AD,temp2)
      CALL myMATMUL(temp2,trans4,temp1)
      CALL myADD(term5_AD,temp1,term5_AD)

      CALL myMATMUL(refl_AD,trans4,temp1)
      CALL myADD(trans_AD,temp1,trans_AD)

      CALL myMATMUL(trans1,refl_AD,temp1)
      CALL myTRANSPOSE(RTV%Trans(:,:,L-1,KL),temp3)
      CALL myMATMUL(temp1,temp3,temp2)
      CALL myADD(refl_AD,temp2,refl_AD)

      CALL myTRANSPOSE(RTV%Refl(:,:,L-1,KL),temp3)
      CALL myMATMUL(term5_AD,temp3,temp1)
      CALL myADD(refl_AD,temp1,refl_AD)

      CALL myMATMUL(temp3,term5_AD, temp1)
      CALL myADD(refl_AD,temp1,refl_AD)

    ENDDO

    s = RTV%Delta_Tau(KL) * single_albedo
    c_AD = ZERO
    s_AD = ZERO
    Delta_Tau_AD=ZERO

!$acc loop independent reduction(+:Delta_Tau_AD,s_AD)
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

    Delta_Tau_AD = Delta_Tau_AD + s_AD* single_albedo
    single_albedo_AD = single_albedo_AD+RTV%Delta_Tau(KL) * s_AD
    optical_depth_AD = optical_depth_AD + Delta_Tau_AD/(TWO**RTV%Number_Doubling(KL))

  END SUBROUTINE CRTM_Doubling_layer_AD

  !------------------------------------------------------------------
  ! print_2d_variable
  !
  ! Prints statistics for a 2d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_2d_variable(name, data)

    character(len=*) :: name
    real(fp)         :: data(:,:)

    ! Note: Assumed shape array sections always start with index=1 for all
    ! dimensions
    !       So we don't have to know start/end indices here
    WRITE(*,'(A5, A17,5ES20.10)') "TEST ", name, minval(data), maxval(data), data(1,1), &
                            data(size(data,1), size(data,2)),            &
                            sqrt(sum(data**2) / size(data))

  END SUBROUTINE print_2d_variable

  !------------------------------------------------------------------
  ! print_3d_variable
  !
  ! Prints statistics for a 3d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_3d_variable(name, data)

    character(len=*) :: name
    real(fp)         :: data(:,:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    WRITE(*,'(A5,A17,5ES20.10)') "TEST ", name, minval(data), maxval(data), data(1,1,1),  &
                            data(size(data,1), size(data,2), size(data,3)), &
                            sqrt(sum(data**2) / size(data))

  END SUBROUTINE print_3d_variable

  !------------------------------------------------------------------
  ! print_4d_variable
  !
  ! Prints statistics for a 4d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_4d_variable(name, data)

    character(len=*) :: name
    real(fp)         :: data(:,:,:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    WRITE(*,'(A5,A17,5ES20.10)') "TEST ", name, minval(data), maxval(data), data(1,1,1,1),  &
                            data(size(data,1), size(data,2), size(data,3), size(data,4)), &
                            sqrt(sum(data**2) / size(data))

  END SUBROUTINE print_4d_variable


  !------------------------------------------------------------------
  ! print_5d_variable
  !
  ! Prints statistics for a 5d state variable
  !------------------------------------------------------------------
  SUBROUTINE print_5d_variable(name, data)

    character(len=*) :: name
    real(fp)         :: data(:,:,:,:,:)

    ! Note: Assumed shape array sections always start with index=1 for all dimensions
    !       So we do not have to know start/end indices here
    WRITE(*,'(A5,A17,5ES20.10)') "TEST ", name, minval(data), maxval(data), data(1,1,1,1,1),  &
                            data(size(data,1), size(data,2), size(data,3), size(data,4), size(data,5)), &
                            sqrt(sum(data**2) / size(data))

  END SUBROUTINE print_5d_variable

  !------------------------------------------------------------------
  ! print_state
  !
  ! Prints statistics for the kernel state variables
  !------------------------------------------------------------------
  SUBROUTINE print_state(msg, &
                MAX_N_ANGLES, &
                    N_LAYERS, &
         N_PROFILESxCHANNELS, &
                      Pff_AD, &
                      Pbb_AD, &
                   s_Refl_AD, &
                  s_Trans_AD, &
              s_source_UP_AD, &
            s_source_DOWN_AD, &
                           w, &
                        T_OD, &
                        w_AD, &
                     T_OD_AD, &
        Planck_Atmosphere_AD, &
                         RTV)

    CHARACTER(LEN=*) :: msg
    INTEGER,        INTENT(IN) :: MAX_N_ANGLES, N_LAYERS, N_PROFILESxCHANNELS
    REAL(fp),       INTENT(IN) :: Pff_AD(:,:,:,:)
    REAL(fp),       INTENT(IN) :: Pbb_AD(:,:,:,:)
    REAL(fp),       INTENT(IN) :: s_Refl_AD(:,:,:,:)
    REAL(fp),       INTENT(IN) :: s_Trans_AD(:,:,:,:)
    REAL(fp),       INTENT(IN) :: s_source_UP_AD(:,:,:)
    REAL(fp),       INTENT(IN) :: s_source_DOWN_AD(:,:,:)
    REAL(fp),       INTENT(IN) :: w(:,:)
    REAL(fp),       INTENT(IN) :: T_OD(:,:)
    REAL(fp),       INTENT(IN) :: w_AD(:,:)
    REAL(fp),       INTENT(IN) :: T_OD_AD(:,:)
    REAL(fp),       INTENT(IN) :: Planck_Atmosphere_AD(:,:)
    TYPE(RTV_type), INTENT(IN) :: RTV(:)

    INTEGER               :: i
    REAL(fp), ALLOCATABLE :: temp2d(:,:), temp3d(:,:,:), temp4d(:,:,:,:), temp5d(:,:,:,:,:)

    WRITE(*,'(A4)') "TEST"
    WRITE(*,'(A5,A117)') "TEST ", repeat("=",117)
    WRITE(*,'(A5,A32)') "TEST ", msg
    WRITE(*,'(A5,A117)') "TEST ", repeat("=",117)
    WRITE(*,'(A5,A17,5A20)') "TEST ", "Variable", "Min", "Max", "First", "Last", "RMS"
    WRITE(*,'(A5,A117)') "TEST ", repeat("-",117)

    CALL print_4d_variable("Pff_AD", Pff_AD)
    CALL print_4d_variable("Pbb_AD", Pbb_AD)
    CALL print_4d_variable("s_Refl_AD", s_Refl_AD)
    CALL print_4d_variable("s_Trans_AD", s_Trans_AD)
    CALL print_3d_variable("s_source_UP_AD", s_source_UP_AD)
    CALL print_3d_variable("s_source_DOWN_AD", s_source_DOWN_AD)
    CALL print_2d_variable("w", w)
    CALL print_2d_variable("T_OD", T_OD)
    CALL print_2d_variable("w_AD", w_AD)
    CALL print_2d_variable("T_OD_AD", T_OD_AD)
    CALL print_2d_variable("Planck_Atmosphere_AD", Planck_Atmosphere_AD)

    ALLOCATE(temp4d(N_PROFILESxCHANNELS, MAX_N_ANGLES, MAX_N_ANGLES + 1, N_LAYERS))
    DO i = 1, N_PROFILESxCHANNELS
      temp4d(i, :, :, :) = RTV(i)%Pff
    END DO
    CALL print_4d_variable("RTV%Pff", temp4d)
    DO i = 1, N_PROFILESxCHANNELS
      temp4d(i, :, :, :) = RTV(i)%Pbb
    END DO
    CALL print_4d_variable("RTV%Pbb", temp4d)
    DEALLOCATE(temp4d)

    ALLOCATE(temp2d(N_PROFILESxCHANNELS, N_LAYERS))
    DO i = 1, N_PROFILESxCHANNELS
      temp2d(i, :) = RTV(i)%Delta_Tau
    END DO
    CALL print_2d_variable("RTV%Delta_Tau", temp2d)
    DEALLOCATE(temp2d)

    ALLOCATE(temp5d(N_PROFILESxCHANNELS, MAX_N_ANGLES, MAX_N_ANGLES, MAX_N_DOUBLING + 1, N_LAYERS))
    DO i = 1, N_PROFILESxCHANNELS
      temp5d(i, :, :, :, :) = RTV(i)%Refl
    END DO
    CALL print_5d_variable("RTV%Refl", temp5d)
    DO i = 1, N_PROFILESxCHANNELS
      temp5d(i, :, :, :, :) = RTV(i)%Trans
    END DO
    CALL print_5d_variable("RTV%Trans", temp5d)
    DO i = 1, N_PROFILESxCHANNELS
      temp5d(i, :, :, :, :) = RTV(i)%Inv_BeT
    END DO
    CALL print_5d_variable("RTV%Inv_BeT", temp5d)
    DEALLOCATE(temp5d)

    ALLOCATE(temp2d(N_PROFILESxCHANNELS, MAX_N_ANGLES))
    DO i = 1, N_PROFILESxCHANNELS
      temp2d(i, :) = RTV(i)%COS_Angle
    END DO
    CALL print_2d_variable("RTV%COS_Angle", temp2d)
    DO i = 1, N_PROFILESxCHANNELS
      temp2d(i, :) = RTV(i)%COS_Weight
    END DO
    CALL print_2d_variable("RTV%COS_Weight", temp2d)
    DEALLOCATE(temp2d)

    ALLOCATE(temp3d(N_PROFILESxCHANNELS, MAX_N_ANGLES, N_LAYERS))
    DO i = 1, N_PROFILESxCHANNELS
      temp3d(i, :, :) = RTV(i)%C1
    END DO
    CALL print_3d_variable("RTV%C1", temp3d)
    DO i = 1, N_PROFILESxCHANNELS
      temp3d(i, :, :) = RTV(i)%C2
    END DO
    CALL print_3d_variable("RTV%C2", temp3d)
    DEALLOCATE(temp3d)

    WRITE(*,'(A5,A117)') "TEST ", repeat("-",117)
    WRITE(*,'(A4)') "TEST"

  END SUBROUTINE print_state

END MODULE my_kernels


!---------------------
!   Driver
!---------------------
PROGRAM test_kernels

  USE my_kernels
  USE omp_lib
#ifdef _OPENACC
  USE openacc
#endif

  !------- Test -------!
  INTEGER, PARAMETER :: N_LAYERS = MAX_N_LAYERS
  INTEGER, PARAMETER :: N_PROFILESxCHANNELS = 100
  TYPE(RTV_type), ALLOCATABLE, DIMENSION(:) :: RTV
  INTEGER :: k, t, alloc_stat, n_omp_threads, streamid
  INTEGER :: N_PROFS_PER_GPU, s, e, gpuid
  INTEGER :: count_rate, count_start, count_end
  REAL :: elapsed

  !---- local arrays --- !
  REAL(fp), ALLOCATABLE, DIMENSION( :,:,:,: ) :: Pff_AD,Pbb_AD,s_Refl_AD,s_Trans_AD
  REAL(fp), ALLOCATABLE, DIMENSION( :,:,: ) :: s_source_UP_AD,s_source_DOWN_AD
  REAL(fp), ALLOCATABLE, DIMENSION(:,:) ::  w, T_OD
  REAL(fp), ALLOCATABLE, DIMENSION(:,:) :: w_AD, T_OD_AD
  REAL(fp), ALLOCATABLE, DIMENSION(:,:) :: Planck_Atmosphere_AD 

  REAL(fp), DIMENSION(MAX_N_ANGLES,MAX_N_ANGLES) :: term1,term2,term3,term4,term5_AD
  REAL(fp), DIMENSION(MAX_N_ANGLES,MAX_N_ANGLES) :: trans1,trans3,trans4,temp1,temp2,temp3
  REAL(fp), DIMENSION(MAX_N_ANGLES) :: C1_AD, C2_AD

  !---- openmp -----!
!$OMP PARALLEL
!$OMP SINGLE
  n_omp_threads = OMP_GET_NUM_THREADS()
!$OMP END SINGLE
!$OMP END PARALLEL

  WRITE(6,*)
  WRITE(6,'("   Using",i3," OpenMP threads for ",i3," profiles and channels.")') &
         n_omp_threads, N_PROFILESxCHANNELS
  WRITE(6,'("   N_LAYERS = ",i3,", N_ANGLES =",i3)') &
         N_LAYERS, MAX_N_ANGLES
  WRITE(6,*)

  !----- multi-gpu ---!
#ifdef _OPENACC
  N_GPUS = acc_get_num_devices(acc_device_nvidia)
#else
  N_GPUS = 0
#endif
  WRITE(6,'(" Number of GPUS = ",i3)') N_GPUS
  IF (N_GPUS .eq. 0) THEN
     N_PROFS_PER_GPU = N_PROFILESxCHANNELS
  ELSE
     N_PROFS_PER_GPU = (N_PROFILESxCHANNELS + N_GPUS - 1) / N_GPUS
  ENDIF

  !---- allocate ----!
  ALLOCATE(RTV(N_PROFILESxCHANNELS), &
           STAT = alloc_stat)
  IF ( alloc_stat /= 0 ) STOP

  ALLOCATE(Pff_AD(MAX_N_ANGLES, MAX_N_ANGLES+1, N_LAYERS, N_PROFILESxCHANNELS), &
           Pbb_AD(MAX_N_ANGLES, MAX_N_ANGLES+1, N_LAYERS, N_PROFILESxCHANNELS), &
           s_Refl_AD(MAX_N_ANGLES, MAX_N_ANGLES, N_LAYERS, N_PROFILESxCHANNELS), &
           s_Trans_AD(MAX_N_ANGLES, MAX_N_ANGLES, N_LAYERS, N_PROFILESxCHANNELS), &
           s_source_UP_AD(MAX_N_ANGLES, N_LAYERS, N_PROFILESxCHANNELS), &
           s_source_DOWN_AD(MAX_N_ANGLES, N_LAYERS, N_PROFILESxCHANNELS), &
           w(N_LAYERS, N_PROFILESxCHANNELS), &
           T_OD(N_LAYERS, N_PROFILESxCHANNELS), &
           w_AD(N_LAYERS, N_PROFILESxCHANNELS), &
           T_OD_AD(N_LAYERS, N_PROFILESxCHANNELS), &
           Planck_Atmosphere_AD(0:N_LAYERS, N_PROFILESxCHANNELS), &
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

  !---- fill RTV ----!
  DO t = 1, N_PROFILESxCHANNELS
      RTV(t)%n_Streams = 6
  END DO

  !--- Copy data to GPU ---!
#ifdef _OPENACC
  DO gpuid = 0, N_GPUS - 1
     CALL acc_set_device_num(gpuid,acc_device_nvidia)
     s = gpuid * N_PROFS_PER_GPU + 1
     e = MIN(N_PROFILESxCHANNELS, s + N_PROFS_PER_GPU - 1)

     WRITE(6,'( "GPU=",i2," section: ",i4,":",i4)') gpuid, s, e

!$acc enter data copyin(Pff_AD(:,:,:,s:e),Pbb_AD(:,:,:,s:e),s_Refl_AD(:,:,:,s:e),s_Trans_AD(:,:,:,s:e), &
!$acc                  s_source_UP_AD(:,:,s:e),s_source_DOWN_AD(:,:,s:e), &
!$acc                  w(:,s:e), T_OD(:,s:e), w_AD(:,s:e), T_OD_AD(:,s:e), Planck_Atmosphere_AD(:,s:e))

  ENDDO
#endif

!$acc wait

  !---- create RTV array ----!
  !---- manual deepcopy with explicit attach ---- !
  !---- For some reason implicit attach does not work ----!
  PRINT*, "Creating RTV"
  DO t = 1, N_PROFILESxCHANNELS
      gpuid = (t - 1) / N_PROFS_PER_GPU
#ifdef _OPENACC
      CALL acc_set_device_num(gpuid,acc_device_nvidia)
#endif
      CALL RTV_Create( RTV(t), MAX_N_ANGLES, MAX_N_LEGENDRE_TERMS, N_LAYERS )
  ENDDO

#ifdef _OPENACC
  DO gpuid = 0, N_GPUS - 1
     CALL acc_set_device_num(gpuid,acc_device_nvidia)
!$acc enter data copyin(RTV)
  ENDDO
#endif

  DO t = 1, N_PROFILESxCHANNELS
      gpuid = (t - 1) / N_PROFS_PER_GPU
#ifdef _OPENACC
      CALL acc_set_device_num(gpuid,acc_device_nvidia)
#endif
!$acc enter data attach(RTV(t)%Pff,RTV(t)%Pbb,RTV(t)%Number_Doubling,RTV(t)%Delta_Tau, &
!$acc                RTV(t)%Refl,RTV(t)%Trans,RTV(t)%Inv_BeT,RTV(t)%C1,RTV(t)%C2)
  ENDDO

  PRINT*, "Finished creating RTV"

  !------- Print input state statistics -------!
  CALL print_state("Input state", &
                    MAX_N_ANGLES, &
                        N_LAYERS, &
             N_PROFILESxCHANNELS, &
                          Pff_AD, &
                          Pbb_AD, &
                       s_Refl_AD, &
                      s_Trans_AD, &
                  s_source_UP_AD, &
                s_source_DOWN_AD, &
                               w, &
                            T_OD, &
                            w_AD, &
                         T_OD_AD, &
            Planck_Atmosphere_AD, &
                             RTV)
  !---- call kernel ----!
  PRINT*, "Calling kernel"  
  CALL SYSTEM_CLOCK (count_rate=count_rate)
  CALL SYSTEM_CLOCK (count=count_start)

!$omp parallel do private(t,k,streamid,gpuid, &
!$omp            term1,term2,term3,term4,term5_AD, &
!$omp            trans1,trans3,trans4,temp1,temp2,temp3,C1_AD,C2_AD)
  DO t = 1, N_PROFILESxCHANNELS

  streamid = mod(t - 1,n_omp_threads)
  gpuid = (t - 1) / N_PROFS_PER_GPU

#ifdef _OPENACC
  CALL acc_set_device_num(gpuid,acc_device_nvidia)
#endif

!$acc kernels async(streamid) &
!$acc present(Pff_AD(:,:,:,t),Pbb_AD(:,:,:,t),s_Refl_AD(:,:,:,t),s_Trans_AD(:,:,:,t), &
!$acc                  s_source_UP_AD(:,:,t),s_source_DOWN_AD(:,:,t), &
!$acc                  w(:,t), T_OD(:,t), w_AD(:,t), T_OD_AD(:,t), Planck_Atmosphere_AD(:,t))
!$acc loop private(k, &
!$acc            term1,term2,term3,term4,term5_AD, &
!$acc            trans1,trans3,trans4,temp1,temp2,temp3,C1_AD,C2_AD)
  DO k = 1, N_LAYERS
       CALL CRTM_Doubling_layer_AD(RTV(t)%n_Streams, RTV(t)%n_Angles, k, w( k, t ), T_OD( k, t ),      &        !Input
                           RTV(t)%COS_Angle, RTV(t)%COS_Weight, RTV(t)%Pff( :, :, k ), RTV(t)%Pbb( :, :, k ), & ! Input
                           RTV(t)%Planck_Atmosphere( k ),    & !Input
                           s_trans_AD( :, :, k, t ), s_refl_AD( :, :, k, t ), s_source_up_AD( :, k, t ),   & 
                           s_source_down_AD( :, k, t ), RTV(t), w_AD( k, t ), T_OD_AD( k, t ), Pff_AD( :, :, k, t ), & 
                           Pbb_AD( :, :, k, t ), Planck_Atmosphere_AD( k, t ), &
                           term1,term2,term3,term4,term5_AD,trans1,trans3,trans4,temp1,temp2,temp3,C1_AD,C2_AD)  !Output
  ENDDO
!$acc end kernels

  ENDDO
!$omp end parallel do

#ifdef _OPENACC
  DO gpuid = 0, N_GPUS - 1
     CALL acc_set_device_num(gpuid,acc_device_nvidia)
!$acc wait
  ENDDO
#endif

  CALL SYSTEM_CLOCK (count=count_end)
  elapsed = REAL (count_end - count_start) / REAL (count_rate)
  PRINT*
  PRINT*
  PRINT*, "Finished executing kernel in =", elapsed  
  PRINT*

  !------- Print section of output -------!
#ifdef _OPENACC
  DO gpuid = 0, N_GPUS - 1
     CALL acc_set_device_num(gpuid,acc_device_nvidia)
     s = gpuid * N_PROFS_PER_GPU + 1
     e = MIN(N_PROFILESxCHANNELS, s + N_PROFS_PER_GPU - 1)
!$acc update self(s_trans_AD(:,:,:,s:e))
!$acc update self(s_refl_AD(:,:,:,s:e))
!$acc update self(s_source_up_AD(:,:,s:e))
!$acc update self(s_source_down_AD(:,:,s:e))
!$acc update self(w(:,s:e))
!$acc update self(T_OD(:,s:e))
!$acc update self(w_AD(:,s:e))
!$acc update self(T_OD_AD(:,s:e))
!$acc update self(Pff_AD(:,:,:,s:e))
!$acc update self(Pbb_AD(:,:,:,s:e))
!$acc update self(Planck_Atmosphere_AD(:,s:e))

  ENDDO
#endif

  PRINT*, s_trans_AD(:,1,1,1)
  !------- Print output state statistics -------!
  CALL print_state("Output state", &
                     MAX_N_ANGLES, &
                         N_LAYERS, &
              N_PROFILESxCHANNELS, &
                           Pff_AD, &
                           Pbb_AD, &
                        s_Refl_AD, &
                       s_Trans_AD, &
                   s_source_UP_AD, &
                 s_source_DOWN_AD, &
                                w, &
                             T_OD, &
                             w_AD, &
                          T_OD_AD, &
             Planck_Atmosphere_AD, &
                              RTV)

END PROGRAM test_kernels
