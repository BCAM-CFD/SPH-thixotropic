!------------------------------------------------------------
!  SUBROUTINE : open files to prepare for writing
!------------------------------------------------------------
      SUBROUTINE io_open(this,rank,&
           num_statis,num_boundary,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_open
        !----------------------------------------------------
        !
        ! Purpose     : Opening files for writting.
        !
        ! Reference   :
        !
        ! Remark      : statistic.dat can be only written
        !               by rank 0.
        !
        ! Revisions   : V0.1 07.12 2009, original version.
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
<
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments :
        !----------------------------------------------------
        
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)		:: rank
        INTEGER,INTENT(IN)              :: num_statis
        INTEGER,INTENT(IN)              :: num_boundary
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_colloid

        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        IF( rank /= 0 ) THEN
           PRINT *,"io_open : ",&
                "File can only be opened by root process !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Open statistics file.
        !----------------------------------------------------
        
        IF ( this%write_output > 0 .AND. &
             num_statis > 0 ) THEN
           
           CALL io_open_statistic(this,rank,stat_info_sub)
           
           IF ( stat_info_sub /= 0) THEN
              PRINT *, "io_open : ", &
                   "opening statistic file has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
#ifdef __IO_COLLOID_SEPARATE
        
        !----------------------------------------------------
        ! Open colloid file for writting individual files.
        !----------------------------------------------------
        
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF( this%write_output > 0 .AND. &
             num_colloid > 0 ) THEN
           
           CAll io_open_colloid(this, &
                rank,num_colloid,stat_info_sub) 
           
           IF ( stat_info_sub /= 0) THEN
              PRINT *, "io_open : ", &
                   "opening colloid file has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
#endif
        
        !----------------------------------------------------
        ! Open boundary file.           
        !----------------------------------------------------
        
        IF ( this%write_output > 0 .AND. &
             num_boundary > 0 ) THEN
           
           CAll io_open_boundary(this,rank,stat_info_sub)             
           
           IF ( stat_info_sub /= 0) THEN
              PRINT *, "io_open : ", &
                   "opening boundary file has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE io_open
      

      SUBROUTINE io_open_statistic(this,rank,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_open_statistic
        !----------------------------------------------------
        !
        ! Purpose     : Opening statistic file for preparing
        !               writting.
        !
        ! Reference   :
        !
        ! Remark      : statistic.dat can be only written
        !               by rank 0.
        !
        ! Revisions   : V0.1 15.07 2009, original version.
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
        ! Arguments.
        !----------------------------------------------------
        
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)		:: rank
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        LOGICAL                         :: read_external 
        LOGICAL                         :: lexist
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        read_external  = &
             control_get_read_external(this%ctrl,stat_info_sub)
        
        IF( rank /= 0) THEN
           PRINT *,"io_open_statistic : "
           PRINT *,"statistic file can only be opened by root process!"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        INQUIRE(FILE=this%statistic_file,EXIST=lexist)
        IF ( lexist ) THEN
           IF( read_external) THEN
              
              OPEN(UNIT=this%statistic_unit,FILE=this%statistic_file,&
                   STATUS="OLD",FORM=this%statistic_fmt,ACTION="WRITE", &
                   POSITION="APPEND",IOSTAT=stat_info_sub)
              PRINT *,"Exisiting statistic file ", &
                   TRIM(this%statistic_file), &
                   " being appended !"  
           ELSE
              PRINT *,"io_open_statistic : "
              PRINT *,"old statistic file ", &
                   TRIM(this%statistic_file), &
                   " is still there !"
              stat_info = -1
              GOTO 9999
           END IF
        ELSE
           OPEN(UNIT=this%statistic_unit,FILE=this%statistic_file,&
                STATUS="NEW",FORM=this%statistic_fmt,ACTION="WRITE", &
                IOSTAT=stat_info_sub)
           PRINT *,"New statistic file ", &
                TRIM(this%statistic_file), " being created !"
        END IF
        
        IF( stat_info_sub /= 0) THEN
           PRINT *, "io_open_statistic :"
           PRINT *, "Creating/Opening the ", &
                TRIM(this%statistic_file), " has error!"
           stat_info = -1
           GOTO 9999
        END IF
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE io_open_statistic

      
      SUBROUTINE io_open_statistic_relax(this,rank,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_open_statistic_relax
        !----------------------------------------------------
        !
        ! Purpose     : Opening statistic file for preparing
        !               writting relax statistic.
        !
        ! Reference   :
        !
        ! Remark      : statistic.dat can be only written
        !               by rank 0.
        !
        ! Revisions   : V0.1 09.03 2010, original version.
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
        ! Arguments.
        !----------------------------------------------------
        
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)		:: rank
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        LOGICAL                         :: read_external 
        LOGICAL                         :: lexist
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        
        read_external  = &
             control_get_read_external(this%ctrl,stat_info_sub)
        
        IF( rank /= 0) THEN
           PRINT *,"io_open_statistic_relax: "
           PRINT *, "statistic_relax file can only be opened by root process!"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        INQUIRE(FILE=this%statistic_relax_file,EXIST=lexist)
        
        IF ( lexist ) THEN
           
           IF( read_external) THEN
              
              OPEN(UNIT=this%statistic_relax_unit, &
                   FILE=this%statistic_relax_file, &
                   STATUS="OLD", &
                   FORM=this%statistic_relax_fmt, &
                   ACTION="WRITE", &
                   POSITION="APPEND",IOSTAT=stat_info_sub)
              
              PRINT *,"Exisiting statistic file ", &
                   TRIM(this%statistic_relax_file), &
                   " being appended !" 
              
           ELSE
              PRINT *,"io_open_statistic_relax : "
              PRINT *,"Old statistic_relax file ", &
                   TRIM(this%statistic_relax_file), &
                   " is still there!"
              stat_info = -1
              GOTO 9999
           END IF
           
        ELSE
           
           OPEN(UNIT=this%statistic_relax_unit, &
                FILE=this%statistic_relax_file,&
                STATUS="NEW", &
                FORM=this%statistic_relax_fmt, &
                ACTION="WRITE", &
                IOSTAT=stat_info_sub)
           
           PRINT *,"New statistic_relax file ", &
                TRIM(this%statistic_relax_file), &
                " being created !"
           
        END IF
        
        IF( stat_info_sub /= 0) THEN
           PRINT *,"io_open_statistic_relax :"
           PRINT *,"creating/opening the ", &
                TRIM(this%statistic_relax_file), " has error!"
           stat_info = -1
           GOTO 9999
        END IF
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE io_open_statistic_relax
      
      
      SUBROUTINE io_open_colloid(this,rank,num_colloid,stat_info)
        !----------------------------------------------------
        !  Subroutine   :  io_open_colloid
        !----------------------------------------------------
        !
        !  Purpose      :  Opening colloid file for preparing
        !                  writting.
        !
        !  Reference    :
        !
        !  Remark       : colloid*.dat can be only written
        !                by rank 0, 1, or 2, only one of
        !                them, not multiple.
        !
        !  Revisions    : V0.1 15.07.2009, original version.
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
        
        !--------------------------------
        ! Arguments.
        !--------------------------------
        
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)		:: rank
        INTEGER, INTENT(IN)             :: num_colloid
        INTEGER, INTENT(OUT)            :: stat_info
        
        !--------------------------------
        ! Local variables.
        !--------------------------------
        
        LOGICAL                         :: read_external
        LOGICAL                         :: lexist
        INTEGER                         :: i
        CHARACTER(len=MAX_CHAR)         :: file_name
        INTEGER                         :: stat_info_sub

        !--------------------------------
        ! Initialization of variables.
        !--------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        read_external  = &
             control_get_read_external(this%ctrl,stat_info_sub)
        
        IF( rank /= 0 ) THEN
           
           PRINT *,"io_open_colloid : "
           PRINT *,"colloid file can only be opened by rank 0 !"
           stat_info = -1
           GOTO 9999
           
        ELSE
           
           DO i = 1, num_colloid
              
              WRITE(file_name,'(A,I4.4,A)') & 
                   TRIM(this%colloid_file),i,'.dat'                   
              
              INQUIRE(FILE=file_name,EXIST=lexist)
              
              IF ( lexist ) THEN
                 
                 IF( read_external) THEN
                    PRINT *,"Exisiting colloid file ", &
                         TRIM(file_name), &
                         " being appended !" 
                    
                    OPEN(UNIT=this%colloid_unit+i,FILE=file_name,&
                         STATUS="OLD",FORM=this%colloid_fmt,ACTION="WRITE", &
                         POSITION="APPEND",IOSTAT=stat_info_sub)
                 ELSE
                    
                    PRINT *,"io_open_colloid   : "
                    PRINT *,"Old colloid file ",&
                         TRIM(file_name), &
                         "is still there !"
                    
                    stat_info = -1
                    GOTO 9999
                    
                 END IF
              ELSE
                 PRINt *,"New colloid file ", &
                      TRIM(file_name),&
                      " being created !"
                 
                 OPEN(UNIT=this%colloid_unit+i, &
                      FILE=file_name, &
                      STATUS="NEW", &
                      FORM=this%colloid_fmt, &
                      ACTION="WRITE", &
                      IOSTAT=stat_info_sub)
              END IF
              
              IF(stat_info_sub /= 0) THEN
                 PRINT *,"io_open_colloid :"
                 PRINT *,"Creating/Opening file ", &
                      TRIM(file_name), " has error!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END DO
           
        END IF
           
           
9999    CONTINUE	
        
        
        RETURN
        
      END SUBROUTINE io_open_colloid 
      

      SUBROUTINE io_open_boundary(this,rank,stat_info)
        !----------------------------------------------------
        !  Subroutine   :  io_open_boundary
        !----------------------------------------------------
        !
        !  Purpose      :  Opening boundary file for preparing
        !                  writting.
        !
        !  Reference    :
        !
        !  Remark       :  boundary*.dat can be only written
        !                  by rank 0 or 1, only one of
        !                  them, not multiple.
        !
        !  Revisions    : V0.1 15.07.2009, original version.
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
        
        !--------------------------------
        ! Arguments.
        !--------------------------------
        
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)		:: rank
        INTEGER, INTENT(OUT)            :: stat_info
        
        !--------------------------------
        ! Local variables.
        !--------------------------------
        
        LOGICAL                         :: read_external 
        LOGICAL                         :: lexist
        INTEGER                         :: stat_info_sub
        
        !--------------------------------
        ! Initialization of variables.
        !--------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        read_external  = &
             control_get_read_external(this%ctrl,stat_info_sub)
        
        IF( rank > 1 ) THEN
           PRINT *,"io_open_bounary : "
           PRINT *,"boundary file can only be opened by rank 0, or 1 !"
           stat_info = -1
           GOTO 9999           
        ELSE
           
           INQUIRE(FILE=this%boundary_file,EXIST=lexist)
           IF ( lexist ) THEN
              IF( read_external) THEN
                 
                 OPEN(UNIT=this%boundary_unit,FILE=this%boundary_file,&
                      STATUS="OLD",FORM=this%boundary_fmt,ACTION="WRITE", &
                      POSITION="APPEND",IOSTAT=stat_info_sub)
                 PRINT *,"Exisiting boundary file ",&
                      TRIM(this%boundary_file), &
                      " being appended !"  
              ELSE
                 PRINT *,"io_open_boundary  : "
                 PRINT *, "Old boundary file ",&
                      TRIM(this%boundary_file), &
                      " is still there !"
                 stat_info = -1
                 GOTO 9999
              END IF
           ELSE
              OPEN(UNIT=this%boundary_unit,FILE=this%boundary_file,&
                   STATUS="NEW",FORM=this%boundary_fmt,ACTION="WRITE", &
                   IOSTAT=stat_info_sub)
              PRINT *,"New boundary file ",&
                      TRIM(this%boundary_file), &
                      " being created !"
           END IF
           
           IF( stat_info_sub /= 0) THEN
              PRINT *, "io_open_boundary :"
              PRINT *, "creating/opening the ", &
                   TRIM(this%boundary_file), " has error!"              
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE io_open_boundary
      
