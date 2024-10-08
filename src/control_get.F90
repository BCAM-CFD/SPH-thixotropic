!--------------------------------------------------
! Subroutine  :  control_get_*
!--------------------------------------------------
!
! Purpose     : Get routines of Class control.
!
! Reference   :
!
! Remark      :
!
! Revisions   : V0.1 01.03.2009, original version.
!
!--------------------------------------------------
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
	!--------------------------------------------------

      INTEGER FUNCTION control_get_debug_flag(this,stat_info)

        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        control_get_debug_flag = this%debug_flag
        
        RETURN 
        
      END FUNCTION control_get_debug_flag
      
      
      LOGICAL FUNCTION control_get_relax_run(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_relax_run  = this%relax_run
        
        RETURN
        
      END FUNCTION control_get_relax_run

      
      LOGICAL FUNCTION control_get_colloid_relax(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_colloid_relax  = this%colloid_relax
        
        RETURN
        
      END FUNCTION control_get_colloid_relax
      
      
      LOGICAL FUNCTION control_get_read_external(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_read_external  = this%read_external
        
        RETURN
        
      END FUNCTION control_get_read_external
      
      
      INTEGER FUNCTION control_get_kernel_type(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_kernel_type = this%kernel_type
        
        RETURN
        
      END FUNCTION control_get_kernel_type
      
      
      LOGICAL FUNCTION control_get_symmetry(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_symmetry  = this%symmetry
        
        RETURN
        
      END FUNCTION control_get_symmetry     
      
      
      INTEGER FUNCTION control_get_rhs_density_type(this,stat_info)

        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info

        stat_info = 0
        control_get_rhs_density_type = this%rhs_density_type
        
        RETURN
        
      END FUNCTION control_get_rhs_density_type
    

      LOGICAL FUNCTION control_get_dynamic_density_ref(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_dynamic_density_ref  = this%dynamic_density_ref
        
        RETURN
        
      END FUNCTION control_get_dynamic_density_ref
      

      INTEGER FUNCTION control_get_stateEquation_type(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        control_get_stateEquation_type = this%stateEquation_type
        
        RETURN

      END FUNCTION control_get_stateEquation_type
      

      LOGICAL FUNCTION control_get_Newtonian(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_Newtonian  = this%Newtonian
        
        RETURN
        
      END FUNCTION control_get_Newtonian


      LOGICAL FUNCTION control_get_Brownian(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_Brownian  = this%Brownian
        
        RETURN
        
      END FUNCTION control_get_Brownian


      INTEGER FUNCTION control_get_random_seed(this,stat_info)

        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info

        stat_info = 0
        control_get_random_seed = this%random_seed
        
        RETURN
        
      END FUNCTION control_get_random_seed
      

      INTEGER FUNCTION control_get_rhs_force_type(this,stat_info)

        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info

        stat_info = 0
        control_get_rhs_force_type = this%rhs_force_type
        
        RETURN
        
      END FUNCTION control_get_rhs_force_type

      
      LOGICAL FUNCTION control_get_pp_interact_cc(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_pp_interact_cc  = this%pp_interact_cc
        
        RETURN
        
      END FUNCTION control_get_pp_interact_cc

      
      LOGICAL FUNCTION control_get_pp_interact_cw(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_pp_interact_cw  = this%pp_interact_cw
        
        RETURN
        
      END FUNCTION control_get_pp_interact_cw
      
      
      INTEGER FUNCTION control_get_cc_lub_type(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_cc_lub_type  = this%cc_lub_type
        
        RETURN
        
      END FUNCTION control_get_cc_lub_type


      INTEGER FUNCTION control_get_cc_repul_type(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_cc_repul_type  = this%cc_repul_type
        
        RETURN
        
      END FUNCTION control_get_cc_repul_type


      INTEGER FUNCTION control_get_cc_magnet_type(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_cc_magnet_type  = this%cc_magnet_type
        
        RETURN
        
      END FUNCTION control_get_cc_magnet_type

      
      INTEGER FUNCTION control_get_cw_lub_type(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_cw_lub_type  = this%cw_lub_type
        
        RETURN
        
      END FUNCTION control_get_cw_lub_type
      
      
      INTEGER FUNCTION control_get_cw_repul_type(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_cw_repul_type  = this%cw_repul_type
        
        RETURN
        
      END FUNCTION control_get_cw_repul_type
      
      
      LOGICAL FUNCTION control_get_stress_tensor(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_stress_tensor  = this%stress_tensor
        
        RETURN
        
      END FUNCTION control_get_stress_tensor

      
      LOGICAL FUNCTION control_get_p_energy(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_p_energy  = this%p_energy
        
        RETURN
        
      END FUNCTION control_get_p_energy

      
      LOGICAL FUNCTION control_get_flow_v_fixed(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_flow_v_fixed  = this%flow_v_fixed
        
        RETURN
        
      END FUNCTION control_get_flow_v_fixed
      
      
      INTEGER FUNCTION control_get_integrate_type(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
       
        stat_info = 0
        
        control_get_integrate_type = this%integrate_type
       
        RETURN
       
      END FUNCTION control_get_integrate_type

      
      INTEGER FUNCTION control_get_integrate_colloid_type(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
       
        stat_info = 0
        
        control_get_integrate_colloid_type = this%integrate_colloid_type
       
        RETURN
       
      END FUNCTION control_get_integrate_colloid_type


      INTEGER FUNCTION control_get_integrate_colloid_RK(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
       
        stat_info = 0
        
        control_get_integrate_colloid_RK = this%integrate_colloid_RK
       
        RETURN
       
      END FUNCTION control_get_integrate_colloid_RK


      INTEGER FUNCTION control_get_integrate_colloid_AB(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
       
        stat_info = 0
        
        control_get_integrate_colloid_AB = this%integrate_colloid_AB
       
        RETURN
       
      END FUNCTION control_get_integrate_colloid_AB

      
      INTEGER FUNCTION control_get_adaptive_dt(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_adaptive_dt  = this%adaptive_dt
        
        RETURN
        
      END FUNCTION control_get_adaptive_dt

      
      INTEGER FUNCTION control_get_write_output(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_write_output  = this%write_output
        
        RETURN
        
      END FUNCTION control_get_write_output


      INTEGER FUNCTION control_get_write_restart(this,stat_info)
        
        TYPE(Control), INTENT(IN)                :: this
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        control_get_write_restart  = this%write_restart
        
        RETURN
        
      END FUNCTION control_get_write_restart

      
      REAL(MK) FUNCTION control_get_elapsed_wall_time(this,stat_info)

        !####################################################
        ! return elapsed wall time in hours
        ! since this run starts.
        !####################################################

        TYPE(Control), INTENT(IN)       :: this
        INTEGER, INTENT(out)            :: stat_info
        
        REAL(MK)                        :: elapse
        
        stat_info = 0
        
        CALL CPU_TIME(elapse)
        
        control_get_elapsed_wall_time = elapse / 3600.0_MK
        
        RETURN
        
      END FUNCTION control_get_elapsed_wall_time
