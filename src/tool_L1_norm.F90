      REAL(MK) FUNCTION tool_L1_norm(this,x,stat_info)
        !----------------------------------------------------
        ! Subroutine  : tool_L1_norm
        !----------------------------------------------------
        !
        ! Purpose     : Calculate L1 norm of a vector.
        !
        ! Routines    :
        !
        ! Remarks     : 
        !
        ! References  :
        !
        ! Revisions   : V0.1 25.05 2012, original version.
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
        ! Arguments     
        !----------------------------------------------------
        
        TYPE(Tool), INTENT(IN)                  :: this
        REAL(MK), DIMENSION(:),INTENT(IN)       :: x
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: num, i
        REAL(MK)                                :: norm
        
        !----------------------------------------------------
        ! Initialization
        !
        ! This is supposed to be used, otherwise,
        ! compiler complains that it is not used.
        !----------------------------------------------------
        
        stat_info = this%flag
        stat_info = 0
        
        
        num = SIZE(x,1)
        
        norm = 0.0_MK
        
        DO i = 1, num
           
           norm  = norm + ABS(x(i))
           
        END DO
        
        tool_L1_norm = norm
        
        !----------------------------------------------------
        ! Return.
        !----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END FUNCTION  tool_L1_norm
      
      
