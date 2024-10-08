	!----- Added by Adolfo -------
  !----------------------------------------------------------
  ! particle-particle_interaction_wall_fluid
  !
  ! ip is wall boundary particle. 
  ! jp is fluid particle.
  !
  ! In case of solid wall boundary,
  ! substitute ip's velocity with an artificial
  ! velocity for no slip condition,
  ! pressure ip with Neumann condition.
  !
  ! If Non-Newtonian fluids, 
  ! calculate velocity gradient tensor.
  ! pressure tensor of jp with Neumann condition.
  !
  !
  ! w_f interaction takes place the same way
  ! as f_f interaction.
  !----------------------------------------------------------
  

  v_ip(1:num_dim) = this%v(1:num_dim,ip)
  v_jp(1:num_dim) = this%v(1:num_dim,jp)

  IF ( num_wall_solid > 0 ) THEN
     
     !-------------------------------------------------------
     ! wall boundary using solid particles with
     ! Morris no slip boundary condition.
     !-------------------------------------------------------
   
     IF ( wall_noslip == 2  ) THEN
        
        CALL boundary_noslip(tboundary, &
             this%x(1:num_dim,jp),this%x(1:num_dim,ip), &
             this%v(1:num_dim,jp),v_ip(1:num_dim), &
             this%id(this%sid_idx,ip), &
             stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN
           PRINT *, "pp_interaction_wf : ", &
                "Wall no slip of jp has problem !"     
           stat_info = -1
           GOTO 9999
        END IF
        
     END IF  ! wall_noslip
     
  END IF ! num_wall_solid > 0
  
     !-------------------------------------------------------
     ! For non-Newtonian fluid, density and pressure 
     ! using Neumann conditon.
     !-------------------------------------------------------

     rho_ip = this%rho(jp)
     rho_jp = this%rho(jp)
     
#include "pp_vgt_ip.inc"
     
     IF ( symmetry ) THEN
        
#include "pp_vgt_jp.inc"
        
     END IF
