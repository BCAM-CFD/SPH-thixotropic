      SUBROUTINE particles_compute_act_OLD(this,num,stat_info)
        !----------------------------------------------------
        ! Subroutine :  particles_compute_act
        !----------------------------------------------------
        !
        ! Purpose   :  Compute the accleration of conformation
        !              tensor of particles.
        !
        ! Reference : Vazquez-Quesada, Ellero and Espanol
        !             Phyical Review E 79. 056707, 2009.
        !
        ! Remark    :
        !
        ! Revision  :  V0.1  31.07.2009, original version.
        !
        !              15.04.2010, Adolfo Vázquez-Quesada
        !                 I changed the coefficients of the 
        !                 velocity gradient tensor in act
        !----------------------------------------------------
        ! This code is  based on the original MCF code  developed by Xin Bian.
        ! The  current version  has  been developed  in collaboration  between
        ! - Marco Ellero,  leader of the  CFD Modelling and Simulation  group at
        ! BCAM (Basque Center  for Applied Mathematics) in  Bilbao, Spain
	! - Emanuele Rossi, from Marco Ellero's group.
        ! - Adolfo Vazquez-Quesada from  the Department of Fundamental Physics
        ! at UNED, in Madrid, Spain.
        !
        ! Developers:
        !     Xin Bian.
        !     Adolfo Vazquez-Quesada.
	!     Emanuele Rossi
        !
        ! Contact: a.vazquez-quesada@fisfun.uned.es
        !          mellero@bcamath.org

        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !
        ! this           : an object of Particles Class.
        ! num            : number of particles updated,
        !                  i.e. first num particles in this%x 
        !                  are operated.
        ! lambda         : coefficient required.
        ! stat_info      : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        REAL(MK)                                :: tau
        INTEGER                                 :: i,j,k,p
        REAL(MK), DIMENSION(3,3)                :: t_c
        REAL(MK), DIMENSION(3,3)                :: t_vgt
        REAL(MK), DIMENSION(3,3)                :: Im
        REAL(MK)                                :: t_act

        !-------------------------------
        ! Initialization of variables.
        !-------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        !--------------------------------
        ! Calculation only for real particles.
        !--------------------------------
        
        IF( num > this%num_part_real) THEN
           PRINT *, "particles_compute_act : ", &
                "num > num_part_real, wrong !"
           stat_info = -1
           GOTO 9999      
        END IF
        
        dim = physics_get_num_dim(this%phys,stat_info_sub)
        tau = physics_get_tau(this%phys,stat_info_sub)
        
        IF (ABS(tau) < mcf_machine_zero) THEN
           PRINT *, "particles_compute_act : ", &
                "tau should not be zero !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !-------------------------------
        ! Allocate memory for act.
        !-------------------------------
        
        IF(ASSOCIATED(this%act)) THEN
           DEALLOCATE(this%act)
        ENd IF
        
        ALLOCATE(this%act(dim**2,num))
        
        this%act(:,:) = 0.0_MK
        
        !----------------------
        ! Identity matrix
        !----------------------
        
        Im(1:3,1:3) = 0.0_MK
        Im(1,1) = 1.0_MK
        Im(2,2) = 1.0_MK
        Im(3,3) = 1.0_MK
        
        !------------------------------------------
        ! Compute accleration of conformation 
        ! tensor using velocity gradient tensor.
        !------------------------------------------
        
        DO p = 1, num
           
           !---------------------------------------
           ! Convert array notation to matrix
           ! for clarity.
           !---------------------------------------
           
           DO j = 1, dim              
              DO i = 1, dim
                 t_c(i,j)   = this%ct(i+dim*(j-1),p)
                 t_vgt(i,j) = this%vgt(i+dim*(j-1),p)
              END DO
           END DO
                      
           DO j =1, dim  !--- row direction
              
              Do i=1, dim ! | column direction
                 
                 t_act = 0.0_MK
                 
                 DO k= 1, dim
                    
                    !------------------------------
                    ! c <dot> k + k^t <dot> c
                    !------------------------------
                    
!!$                    t_act = t_act + &
!!$                         t_c(i,k) * t_vgt(k,j) + &
!!$                         t_vgt(k,i) * t_c(k,j)

                    t_act = t_act + &
                         t_c(i,k) * t_vgt(j,k) + &
                         t_vgt(i,k) * t_c(k,j)
                    
                 END DO ! k
                 
                 
                 this%act(i+dim*(j-1),p) = t_act - &
                      ( t_c(i,j) - Im(i,j) ) /tau                    
                 
              END DO ! i
              
           END DO ! j
           
        END DO ! p
        
9999    CONTINUE      
        
        RETURN
        
      END SUBROUTINE particles_compute_act_OLD


!----------- Added by Emanuele: thyxotropic model      
      SUBROUTINE particles_compute_act(this,num,stat_info)
        !----------------------------------------------------
        ! Subroutine :  particles_compute_act
        !----------------------------------------------------
        !
        ! Purpose   :  Compute the accleration of conformation
        !              tensor of particles.
        !
        ! Reference : Vazquez-Quesada, Ellero and Espanol
        !             Phyical Review E 79. 056707, 2009.
        !
        ! Remark    :
        !
        ! Revision  :  V0.1  31.07.2009, original version.
        !
        !              15.04.2010, Adolfo Vázquez-Quesada
        !                 I changed the coefficients of the 
        !                 velocity gradient tensor in act
        !----------------------------------------------------
        ! This code is  based on the original MCF code  developed by Xin Bian.
        ! The  current version  has  been developed  in collaboration  between
        ! - Marco Ellero,  leader of the  CFD Modelling and Simulation  group at
        ! BCAM (Basque Center  for Applied Mathematics) in  Bilbao, Spain
	! - Emanuele Rossi, from Marco Ellero's group.
        ! - Adolfo Vazquez-Quesada from  the Department of Fundamental Physics
        ! at UNED, in Madrid, Spain.
        !
        ! Developers:
        !     Xin Bian.
        !     Adolfo Vazquez-Quesada.
	!     Emanuele Rossi
        !
        ! Contact: a.vazquez-quesada@fisfun.uned.es
        !          mellero@bcamath.org

        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !
        ! this           : an object of Particles Class.
        ! num            : number of particles updated,
        !                  i.e. first num particles in this%x 
        !                  are operated.
        ! lambda         : coefficient required.
        ! stat_info      : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        REAL(MK)                                :: tau
        INTEGER                                 :: i,j,k
        REAL(MK), DIMENSION(3,3)                :: t_c
        REAL(MK), DIMENSION(3,3)                :: t_vgt
        REAL(MK), DIMENSION(3,3)                :: Im
        REAL(MK)                                :: t_act
        
        !----------------------------------------------------
        INTEGER :: a,b,p
        
        !-------------------- vgt variables -------------
        REAL(MK) :: gamma_dot
        REAL(MK), DIMENSION(:), ALLOCATABLE :: vgt_sym
        !------------------------------------------------

        !-------------------------------
        ! Initialization of variables.
        !-------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        !--------------------------------
        ! Calculation only for real particles.
        !--------------------------------
        
        IF( num > this%num_part_real) THEN
           PRINT *, "particles_compute_act : ", &
                "num > num_part_real, wrong !"
           stat_info = -1
           GOTO 9999      
        END IF
        
        dim = physics_get_num_dim(this%phys,stat_info_sub)
        tau = physics_get_tau(this%phys,stat_info_sub)
        
        IF (ABS(tau) < mcf_machine_zero) THEN
           PRINT *, "particles_compute_act : ", &
                "tau should not be zero !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !-------------------------------
        ! Allocate memory for act.
        !-------------------------------
        
        IF(ASSOCIATED(this%act)) THEN
           DEALLOCATE(this%act)
        ENd IF
        
        ALLOCATE(this%act(dim**2,num))
        
        !-- vgt_sym ---
        ALLOCATE(vgt_sym(dim**2))
        !-------------- 
        
!----------- Added by Emanuele: thyxotropic model ---------------------
! we will use only the first component 
        this%act(:,:) = 0.0_MK
        
        !----------------------
        ! Identity matrix
        !----------------------
        
        Im(1:3,1:3) = 0.0_MK
        Im(1,1) = 0.025_MK
        Im(2,2) = 8.0_MK*Im(1,1)
!~         Im(3,3) = 1.0_MK     
        
        DO p=1,num
           DO b = 1, dim  ! ---, row direction
              DO a = 1, dim ! |,  column direction                          
                 vgt_sym(a+dim*(b-1)) = &
                      this%vgt(a+dim*(b-1), p) + &
                      this%vgt(b+dim*(a-1), p)
              ENDDO
           ENDDO

           !-- The shear rate is calculated --
           IF (dim == 2) THEN
              gamma_dot = &
                   vgt_sym(1)*vgt_sym(1) + &
                   vgt_sym(4)*vgt_sym(4) + &
                   2.0_MK * vgt_sym(2)*vgt_sym(3) 
           ELSE ! dim = 3
              gamma_dot = &
                   vgt_sym(1)*vgt_sym(1) + &
                   vgt_sym(5)*vgt_sym(5) + &
                   vgt_sym(9)*vgt_sym(9) + &
                   2.0_MK * vgt_sym(2)*vgt_sym(4) + &
                   2.0_MK * vgt_sym(3)*vgt_sym(7) + &
                   2.0_MK * vgt_sym(6)*vgt_sym(8)
           ENDIF
           gamma_dot = sqrt(0.5_MK * ABS(gamma_dot))
           
           this%act(1,p) = Im(1,1) - (Im(1,1)+Im(2,2)*gamma_dot)*this%ct(1,p)
        ENDDO
        
        !------------------------------------------
        ! Compute accleration of conformation 
        ! tensor using velocity gradient tensor.
        !------------------------------------------
        
!~         DO p = 1, num
           
!~            !---------------------------------------
!~            ! Convert array notation to matrix
!~            ! for clarity.
!~            !---------------------------------------
           
!~            DO j = 1, dim              
!~               DO i = 1, dim
!~                  t_c(i,j)   = this%ct(i+dim*(j-1),p)
!~                  t_vgt(i,j) = this%vgt(i+dim*(j-1),p)
!~               END DO
!~            END DO
                      
!~            DO j =1, dim  !--- row direction
              
!~               Do i=1, dim ! | column direction
                 
!~                  t_act = 0.0_MK
                 
!~                  DO k= 1, dim
                    
!~                     !------------------------------
!~                     ! c <dot> k + k^t <dot> c
!~                     !------------------------------
                    
!~ !!$                    t_act = t_act + &
!~ !!$                         t_c(i,k) * t_vgt(k,j) + &
!~ !!$                         t_vgt(k,i) * t_c(k,j)

!~                     t_act = t_act + &
!~                          t_c(i,k) * t_vgt(j,k) + &
!~                          t_vgt(i,k) * t_c(k,j)
                    
!~                  END DO ! k
                 
                 
!~                  this%act(i+dim*(j-1),p) = t_act - &
!~                       ( t_c(i,j) - Im(i,j) ) /tau                    
                 
!~               END DO ! i
              
!~            END DO ! j
           
!~         END DO ! p
        
9999    CONTINUE      
        
        RETURN
        
      END SUBROUTINE particles_compute_act
!----------- End of Added by Emanuele: thyxotropic model      
