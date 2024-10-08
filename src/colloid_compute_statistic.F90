      SUBROUTINE  colloid_compute_statistic(this,stat_info)
        !----------------------------------------------------
        !  Program      :   colloid_compute_statistic
        !----------------------------------------------------
        !
        !  Purpose      :   Compute statistics of colloids.
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !
        !  Revisions    :  
        !  Revisions    : V0.2 08.07.2009, 
        !                 check again the work flow is correct and
        !                 supply with more comments for code.
        !
        !                 V0.1 21.06.2009, original version.
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
        !  Arguments
        !
        !  this       : an object of Marching Class.
        !  stat_info  : return flag of status.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(OUT)        :: this
        INTEGER, INTENT(OUT)              :: stat_info
        
      

        INTEGER           :: i,j, dim, num
        REAL(MK)          :: v2
        
        stat_info =  0
        
        dim = this%num_dim
        num = this%num_colloid
        
        this%k_energy(1:num)  = 0.0_MK
        this%mom(1:dim,1:num) = 0.0_MK
        this%k_energy_tot     = 0.0_MK
        this%mom_tot(1:dim)   = 0.0_MK
        
        DO i = 1, num
           
           v2 = 0.0_MK
           
           DO j=1, dim 
              
              v2= v2+ this%v(j,i,1)**2
              
              this%mom(j,i) = this%mom(j,i) + &
                   this%m(i) * this%v(j,i,1)
              
           END DO
           
           this%k_energy(i) = this%k_energy(i) + &
                0.5_MK* this%m(i)*v2
           
           this%k_energy_tot = this%k_energy_tot + &
                this%k_energy(i)
           
           this%mom_tot(1:dim) = this%mom_tot(1:dim) + &
                this%mom(1:dim,i)
           
        END DO ! i = 1, num
        
        
        RETURN
        
      END SUBROUTINE colloid_compute_statistic
