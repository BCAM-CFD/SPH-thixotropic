      SUBROUTINE colloid_create_boundary_particle_3D_dicolloid(this,&
           dx,p_x,sid,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_create_boundary_particle_3D_dicolloid
        !----------------------------------------------------
        !
        ! Purpose     : create the particles on the surface
        !               of a dicolloid.
        !                              
        !
        ! Reference   :
        !
        ! Remark      : 1)
        !               Velocity of boundary particle will
        !               be set to zero, no matter if the
        !               colloid's velocity is not zero,
        !               since it needs to be zero for
        !               relax runs. 
        !               If non-zero velocity(according to
        !               colloid velocity) is needed,
        !               it will be set again after relax
        !               run.
        !
        !
        ! Revisions   : V0.1  Nov. 30, 2011, 
        !               implemented model 5
        !               i.e.,  psfdrm. 
        !
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
        ! Arguments :
        !
        ! Input
        !
        ! this      : object of colloid class.
        ! dx        : initial distance between particles
        ! p_x       : position.
        ! sid       : species ID.
        ! stat_info : status of this routine.
        !----------------------------------------------------
        
        TYPE(colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: dx
        REAL(MK), DIMENSION(:,:), POINTER       :: p_x
        INTEGER, INTENT(IN)                     :: sid
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------
      
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num_dim
        INTEGER                                 :: num_layer
        REAL(MK)                                :: rad,rad_max
        REAL(MK)                                :: d_phi, d_theta
        INTEGER                                 :: num_phi, num_theta
        REAL(MK)                                :: phi, theta, gamma
        INTEGER                                 :: num_max
        INTEGER                                 :: num, num1, num2
        
        INTEGER                                 :: i,j,k,m
        
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: t_x
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        IF ( this%place == mcf_colloid_place_lattice ) THEN
           PRINT *, "colloid_create_boundary_particles_3D_sphere: ", &
                "Boundary particles located on lattice !"
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! Get parameters
        !
        ! num_layer : number of layers of boundary particles
        !             around cutoff thickness.
        !----------------------------------------------------
        
        num_dim = this%num_dim
        
        IF ( num_dim /= 3 ) THEN
           PRINT *, "colloid_create_boundary_particle_3D_sphere : ", &
                "Dimension should be 3 !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        IF ( this%shape(sid) /= mcf_colloid_shape_dicolloid ) THEN
           PRINT *, "colloid_create_boundary_particle_3D_sphere : ", &
                "shape should be sphere !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Number of layers for boundary particles
        !----------------------------------------------------

        num_layer = CEILING( this%din / dx(3) - 0.5_MK)
        
        !PRINT *, "num_layer: ", num_layer
        
        !----------------------------------------------------
        ! Find biggest radius, i.e., outest layer.
        !----------------------------------------------------
        
        rad_max = this%radius(1,sid)
        
        !----------------------------------------------------
        ! Estimate what would be the biggest number of
        ! boundary particles for one colloid.
        ! Each boundary particle's mass is the same and
        ! equal to rho*dx(1)*dx(2)*dx(3)
        !----------------------------------------------------
        
        d_phi   = dx(1) / rad_max
        num_phi = NINT(mcf_pi / d_phi)+2
        
        d_theta   = dx(2) / rad_max
        num_theta = NINT(mcf_pi / d_theta)+2
        
        num_max =  num_layer * num_phi * num_theta
        
        !----------------------------------------------------
        ! To be safe, we double num_max
        !----------------------------------------------------
       
        num_max = num_max * 2
        
        ALLOCATE(t_x(num_dim,num_max))
        
      
        !----------------------------------------------------
        ! Set temporary total number of boundary particles 
        ! to zero.
        !----------------------------------------------------
        
        num = 0
        
        SELECT CASE ( this%place ) 
           
        CASE ( mcf_colloid_place_psfdrm ) 
                 
           DO m = 1, num_layer
                    
              !----------------------------------------------
              ! The first layer is rad - dx(3)/2 away 
              ! from surface, the second is rad-3*dx(3)/2
              !  ...
              !----------------------------------------------
                    
              rad = this%radius(1,sid) - m*dx(3)+ dx(3)/2.0_MK
                    
              IF ( this%radius(2,sid) < rad ) THEN
                 gamma = ACOS(this%radius(2,sid)/rad)
              ELSE
                 gamma = 0.0_MK
              END IF
              
              !PRINT *, "gamma: ", gamma
              
              d_phi   = dx(1) / rad
              num_phi = NINT( (mcf_pi-gamma) / d_phi)
              d_phi   = (mcf_pi-gamma) / DFLOAT(num_phi)
              
              !----------------------------------------------
              ! Right sphere:
              ! The first particle is at east pole.
              !----------------------------------------------
              
              num  = num + 1
              num1 = num
              t_x(1,num) = rad
              t_x(2,num) = 0.0_MK
              t_x(3,num) = 0.0_MK
              
              phi = 0.0_MK

              DO i = 1, num_phi
                 
                 phi       = phi + d_phi
                 d_theta   = dx(2) / ( rad * DSIN(phi))
                 num_theta = NINT(2.0_MK * mcf_pi / d_theta)
                 d_theta   = 2.0_MK * mcf_pi / DFLOAT(num_theta)
                 
                 theta = 0.0_MK
                 
                 DO k = 1, num_theta
                    
                    num = num + 1
                    t_x(1,num) = rad * DCOS(phi)
                    t_x(2,num) = rad * DSIN(phi) * DSIN(theta)
                    t_x(3,num) = rad * DSIN(phi) * DCOS(theta)
                    theta = theta + d_theta
                    
                 END DO ! num_theta
                 
              END DO ! num_phi
              
              !----------------------------------------------
              ! Rotate the boundary particles according to
              ! its orientation.
              ! Translate the boundary particles according
              ! to its center
              !----------------------------------------------
              
              DO j = num1, num
                 t_x(1,j)   = t_x(1,j) +  this%radius(2,sid)
                 t_x(1:3,j) = MATMUL(&
                      this%acc_matrix(1:3,1:3,sid),&
                      t_x(1:3,j))
              END DO
              
              t_x(1,num1:num) = t_x(1,num1:num) + &
                   this%x(1,sid) 
              t_x(2,num1:num) = t_x(2,num1:num) + &
                   this%x(2,sid) 
              t_x(3,num1:num) = t_x(3,num1:num) + &
                   this%x(3,sid)
              
              
              !----------------------------------------------
              ! Left sphere:
              ! The first particle is at west pole.
              !----------------------------------------------
              num  = num + 1
              num2 = num
              
              t_x(1,num) = -rad
              t_x(2,num) = 0.0_MK
              t_x(3,num) = 0.0_MK
              
              phi = 0.0_MK
              
              DO i = 1, num_phi-1
                 
                 phi       = phi + d_phi
                 d_theta   = dx(2) / ( rad * DSIN(phi))
                 num_theta = NINT(2.0_MK * mcf_pi / d_theta)
                 d_theta   = 2.0_MK * mcf_pi / DFLOAT(num_theta)
                 
                 theta = 0.0_MK
                 
                 DO k = 1, num_theta
                    
                    num = num + 1
                    t_x(1,num) = -rad * DCOS(phi)
                    t_x(2,num) = rad * DSIN(phi) * DSIN(theta)
                    t_x(3,num) = rad * DSIN(phi) * DCOS(theta)
                    theta = theta + d_theta
                    
                 END DO ! num_theta
                 
              END DO ! num_phi
            
              !----------------------------------------------
              ! Rotate the boundary particles according to
              ! its orientation.         
              ! Translate the boundary particles according
              ! to its center
              !----------------------------------------------
              
              DO j = num2, num
                 t_x(1,j)   = t_x(1,j) - this%radius(2,sid)                                  
                 t_x(1:3,j) = MATMUL(&
                      this%acc_matrix(1:3,1:3,sid),&
                      t_x(1:3,j))
              END DO
              
              t_x(1,num2:num) = t_x(1,num2:num) + &
                   this%x(1,sid) 
              t_x(2,num2:num) = t_x(2,num2:num) + &
                   this%x(2,sid) 
              t_x(3,num2:num) = t_x(3,num2:num) + &
                   this%x(3,sid)
              
           END DO ! m = 1, num_layer
           
        CASE DEFAULT
           
           PRINT *, "colloid_create_boundary_particle_3D_sphere: ", &
                "placement not available !"
           stat_info = -1
           GOTO 9999
           
        END SELECT ! place
        
        
        !----------------------------------------------------
        ! Allocate memory for output parameters.
        !----------------------------------------------------
        
        IF ( ASSOCIATED(p_x) ) THEN
           DEALLOCATE(p_x)
        END IF
        
        
        ALLOCATE(p_x(num_dim,num))
        
        !----------------------------------------------------
        ! Copy the result to output parameters.
        !----------------------------------------------------
        
        p_x(1:num_dim,1:num) =  t_x(1:num_dim,1:num)
        
        
9999    CONTINUE
        
        
        RETURN
        
      END SUBROUTINE  colloid_create_boundary_particle_3D_dicolloid
      
      SUBROUTINE colloid_create_boundary_particle_3D_ellipsoid(this,&
           dx,p_x,sid,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_create_boundary_particle_3D_ellipsoid
        !----------------------------------------------------
        !
        ! Purpose     : create the particles on the surface
        !               of an ellipsoid
        !                              
        !
        ! Reference   :
        !
        ! Remark      : 1)
        !               Velocity of boundary particle will
        !               be set to zero, no matter if the
        !               colloid's velocity is not zero,
        !               since it needs to be zero for
        !               relax runs. 
        !               If non-zero velocity(according to
        !               colloid velocity) is needed,
        !               it will be set again after relax
        !               run.
        !
        !
        ! Revisions   : V0.1  Nov. 21, 2011, 
        !               implemented model 5
        !               i.e.,  psfdrm. 
        !
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
        ! Arguments :
        !
        ! Input
        !
        ! this      : object of colloid class.
        ! dx        : initial distance between particles
        ! p_x       : position.
        ! sid       : species ID.
        ! stat_info : status of this routine.
        !----------------------------------------------------
        
        TYPE(colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: dx
        REAL(MK), DIMENSION(:,:), POINTER       :: p_x
        INTEGER, INTENT(IN)                     :: sid
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------
        
        REAL(MK), DIMENSION(:), ALLOCATABLE     :: rad
      
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num_dim
        INTEGER                                 :: num_layer
        REAL(MK)                                :: d_phi, d_theta
        INTEGER                                 :: num_phi
        REAL(MK)                                :: phi
        INTEGER                                 :: num_max
        INTEGER                                 :: num
        
        INTEGER                                 :: i,k,m

        
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: t_x, appo_t_x
        REAL(MK), DIMENSION(:), ALLOCATABLE     :: rad_max, angle, length
        REAL(MK) :: f, h, perimeter, half_perimeter, dl, ds
        REAL(MK) :: s, ratio
        INTEGER :: num_arc, npoints, i_arc, iter, inum, iphi
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        IF ( this%place == mcf_colloid_place_lattice ) THEN
           PRINT *, "colloid_create_boundary_particles_3D_sphere: ", &
                "Boundary particles located on lattice !"
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! Get parameters
        !
        ! num_layer : number of layers of boundary particles
        !             around cutoff thickness.
        !----------------------------------------------------
        
        num_dim = this%num_dim
        
        IF ( num_dim /= 3 ) THEN
           PRINT *, "colloid_create_boundary_particle_3D_sphere : ", &
                "Dimension should be 3 !"
           stat_info = -1
           GOTO 9999
        END IF
        
 
        IF ( this%shape(sid) /= mcf_colloid_shape_ellipsoid ) THEN
           PRINT *, "colloid_create_boundary_particle_3D_sphere : ", &
                "shape should be ellipsoid !"
           stat_info = -1
           GOTO 9999
        END IF
        
        ALLOCATE(rad_max(num_dim))
        
        !----------------------------------------------------
        ! Number of layers for boundary particles
        !----------------------------------------------------

        num_layer = CEILING( this%din / dx(3) - 0.5_MK)
        
        !PRINT *, "num_layer: ", num_layer
        
        !----------------------------------------------------
        ! Find biggest radius, i.e., outest layer.
        !----------------------------------------------------
        
        rad_max(:) = this%radius(:,sid)
     
        !----------------------------------------------------
        ! Estimate what would be the biggest number of
        ! boundary particles for one colloid.
        ! Each boundary particle's mass is the same and
        ! equal to rho*dx(1)*dx(2)*dx(3)
        !----------------------------------------------------
        
        f = rad_max(2)/rad_max(1)
        h = ((1.0_MK-f)/(1.0_MK+f))**2
        perimeter = mcf_pi*(rad_max(1)+rad_max(2))*(1.0_MK & 
                    +(3.0_MK*h)/(10.0_MK+sqrt(4.0_MK-3.0_MK*h)) &
                    +(4.0_MK/mcf_pi-14.0_MK/11.0_MK)*h**12)
        half_perimeter = perimeter*0.5_MK
        
        dl = dx(1)
        
        num_arc = NINT(perimeter/dl)
        num_phi = NINT(2*mcf_pi*rad_max(1)/dl)
        
        num_max = num_arc*num_phi*num_layer
        
        !----------------------------------------------------
        ! To be safe, we double num_max
        !----------------------------------------------------
       
        num_max = num_max * 2
        
        ALLOCATE(appo_t_x(num_dim,num_max))
        ALLOCATE(t_x(num_dim,num_max))
        ALLOCATE(angle(num_max))
        ALLOCATE(length(num_max))
      
        !----------------------------------------------------
        ! Set temporary total number of boundary particles 
        ! to zero.
        !----------------------------------------------------
        
        npoints = 0
        
        !----------------------------------------------------
        ! Allocation of radius variable for ellipsoid creation
        !----------------------------------------------------
        
        ALLOCATE(rad(num_dim))
        
        SELECT CASE ( this%place ) 
           
        CASE ( mcf_colloid_place_psfdrm ) 
                 
           DO m = 1, num_layer
               
              !----------------------------------------------
              ! Initialize the counter over particles on
              ! half ellipse
              !----------------------------------------------
               
               num=0
              
              !----------------------------------------------
              ! The first layer is rad - dx(3)/2 away 
              ! from surface, the second is rad-3*dx(3)/2
              !  ...
              !----------------------------------------------
                    
              rad(:) = this%radius(:,sid) - m*dx(:)+ dx(:)/2.0_MK
                    
              !----------------------------------------------
              ! The first particle is at east pole.
              !----------------------------------------------
                    
              num = num + 1
              
              appo_t_x(1,num) = rad(1)
              appo_t_x(2,num) = 0.0_MK
              appo_t_x(3,num) = 0.0_MK
              
              f = rad(2)/rad(1)
              h = ((1-f)/(1+f))**2
              perimeter = mcf_pi*(rad(1)+rad(2))*(1.0_MK& 
                          +(3.0_MK*h)/(10.0_MK+sqrt(4.0_MK-3.0_MK*h)) &
                          +(4.0_MK/mcf_pi-14.0_MK/11.0_MK)*h**12)
              perimeter = perimeter*0.5_MK
              
              dl=dx(1)
              num_arc = NINT(perimeter/dl)
              dl = perimeter/REAL(num_arc)
              ds = mcf_pi*dl/perimeter
              
              s = 0.0_MK
              angle(num) = s
              DO i_arc=1,num_arc-1
                 num = num + 1
                 s = mcf_pi*REAL(i_arc)*dl/perimeter
                 appo_t_x(1,num) = rad(1) * COS(s)
                 appo_t_x(2,num) = rad(2) * SIN(s)
                 appo_t_x(3,num) = 0.0_MK
                 angle(num) = ds
              ENDDO
              
              !----------------------------------------
              ! The last particle is at west pole
              !----------------------------------------           
              
              num = num + 1
              appo_t_x(1,num)= -rad(1)
              appo_t_x(2,num)= 0.0
              appo_t_x(3,num)= 0.0
              angle(num)=mcf_pi
              
              !----------------------------------------
              ! Iterations to adjust particle positions
              ! The result are equispaced points on
              ! ellipse perimeter
              !----------------------------------------
              
              DO iter=1,10
                 length(1)=0.0
                 DO inum=2,num
                    length(inum)=sqrt((appo_t_x(1,inum)-appo_t_x(1,inum-1))**2+(appo_t_x(2,inum)-appo_t_x(2,inum-1))**2)
                    ratio = (dl-length(inum))/dl
                    angle(inum) = (1.0_MK+ratio)*angle(inum)
                 ENDDO
                 length(num)=0.0
                 
                 s=0.0_MK
                 appo_t_x(1,1)=rad(1)
                 appo_t_x(2,1)=0.0_MK
                 appo_t_X(3,1)=0.0_MK
                 
                 DO inum=2,num-1
                    s=s+angle(inum)
                    appo_t_x(1,inum)=rad(1)*COS(s)
                    appo_t_x(2,inum)=rad(2)*SIN(s)
                    appo_t_x(3,inum)=0.0_MK
                 ENDDO
                 
                 appo_t_x(1,num)=-rad(1)
                 appo_t_x(2,num)=0.0_MK
                 appo_t_X(3,num)=0.0_MK
                 
              ENDDO
              
              !--------------------------------------------
              ! Place equispaced points on ellipsoid
              ! surface
              !--------------------------------------------
              
              npoints = npoints + 1
              
              t_x(1,npoints) = rad(1)
              t_x(2,npoints) = 0.0_MK
              t_x(3,npoints) = 0.0_MK
              
              phi =0.0_MK
              
              DO inum=2,num-1
                 num_phi=NINT(2.0_MK*mcf_pi*appo_t_x(2,inum)/dl) 
                 
                 d_phi = 2.0_MK*mcf_pi/REAL(num_phi)
                 
                 DO iphi=1,num_phi
                    npoints = npoints + 1
                    t_x(1,npoints) = appo_t_x(1,inum)
                    t_x(2,npoints) = appo_t_x(2,inum)*SIN(phi)
                    t_x(3,npoints) = appo_t_x(2,inum)*COS(phi)
                    phi = phi + d_phi
                 ENDDO
              ENDDO
              
              npoints = npoints + 1

              t_x(1,npoints) = -rad(1)
              t_x(2,npoints) = 0.0_MK
              t_x(3,npoints) = 0.0_MK
                 
           END DO ! m = 1, num_layer
           
        CASE DEFAULT
           
           PRINT *, "colloid_create_boundary_particle_3D_sphere: ", &
                "placement not available !"
           stat_info = -1
           GOTO 9999
           
        END SELECT ! place
        
        
        !----------------------------------------------------
        ! Rotate each boundary particle according to initial data
        !----------------------------------------------------
        
        DO i=1,npoints
           t_x(1:3,i) = MATMUL(this%acc_matrix(1:3,1:3,sid),t_x(1:3,i))
        ENDDO
        
        !----------------------------------------------------
        ! Translate each boundary particle with its
        ! colloid center.
        !---------------------------------------------------- 
        
        DO i = 1, num_dim
              
           t_x(i,1:npoints) = t_x(i,1:npoints) + this%x(i,sid)
           
        END DO
        
        
        !----------------------------------------------------
        ! Allocate memory for output parameters.
        !----------------------------------------------------
        
        IF ( ASSOCIATED(p_x) ) THEN
           DEALLOCATE(p_x)
        END IF
        
        
        ALLOCATE(p_x(num_dim,npoints))
        
        !----------------------------------------------------
        ! Copy the result to output parameters.
        !----------------------------------------------------
        
        p_x(1:num_dim,1:npoints) = t_x(1:num_dim,1:npoints)

9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  colloid_create_boundary_particle_3D_ellipsoid

      
