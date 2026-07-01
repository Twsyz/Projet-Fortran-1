! cd /c/Users/moi/Desktop/buro/N7/Methodes\ Volumes\ Finis/BE/ADV\1D
PROGRAM main_total
    USE m_type
    IMPLICIT NONE
    
    TYPE(phys)        :: recup
    TYPE(num)         :: get
    TYPE(taille)      :: dim
    TYPE(phys_time)   :: time

    REAL, DIMENSION(:),   ALLOCATABLE :: x_reg, y_reg
    REAL, DIMENSION(:,:), ALLOCATABLE :: x_noeud, y_noeud
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: U
    REAL, DIMENSION(:,:), ALLOCATABLE :: U_face_x  ! (Nx+1, Ny)
    REAL, DIMENSION(:,:), ALLOCATABLE :: V_face_y  ! (Nx, Ny+1)
    REAL :: V(2)
    REAL, DIMENSION(:,:), ALLOCATABLE :: C, C_new

    REAL :: Cinit,MIN_VALUE

    INTEGER :: i,j,n
    REAL :: t, dt

    ! -----------------------------
    ! Lecture des données
    ! -----------------------------
    CALL read_data(recup, get, dim, time)

    ! -----------------------------
    ! Allocation
    ! -----------------------------
    ALLOCATE(x_reg(get%Nx))
    ALLOCATE(y_reg(get%Ny))
    ALLOCATE(x_noeud(get%Nx+1, get%Ny+1))
    ALLOCATE(y_noeud(get%Nx+1, get%Ny+1))
    ALLOCATE(U(get%Nx, get%Ny, 2))
    ALLOCATE(U_face_x(get%Nx+1, get%Ny))
    ALLOCATE(V_face_y(get%Nx, get%Ny+1))
    ALLOCATE(C(get%Nx,get%Ny))
    ALLOCATE(C_new(get%Nx,get%Ny))

    PRINT *, 'Allocations réussies'
    PRINT *, 'get%Nx =', get%Nx, 'get%Ny =', get%Ny
    


    ! -----------------------------
    ! Construction du maillage
    ! -----------------------------
    CALL mesh(x_reg, y_reg, dim, get)
    CALL noeud(x_noeud, y_noeud, dim, get) 
    ! -----------------------------
    ! Initialisation de la concentration
    ! -----------------------------
    DO j = 1, get%Ny
        DO i = 1, get%Nx
            C(i,j) = Cinit(recup, get, dim, x_reg(i), y_reg(j))
        END DO
    END DO
   
    ! -----------------------------
    ! Calcul du pas de temps
    ! -----------------------------

    CALL Time_Step(recup, dim, get, time, x_reg, y_reg, dt)
    
    ! -----------------------------
    ! Calcul final
    ! -----------------------------
    ! OPEN(20, FILE="C_all.dat", STATUS="REPLACE")

    ! WRITE(20,*) get%Nx, get%Ny, get%Nt
    ! WRITE(20,*) dim%L, dim%D, time%tf
    ! WRITE(20,*) recup%alpha, recup%K, time%CFL, time%R

    t = 0.0
    n = 0   
    CALL Champ_Vitesse_Faces(recup, dim, get, x_noeud, y_noeud, U_face_x, V_face_y)

    PRINT *, 'Dimensions de U_face_x :', LBOUND(U_face_x,1), UBOUND(U_face_x,1), LBOUND(U_face_x,2), UBOUND(U_face_x,2)
    PRINT *, 'Dimensions de V_face_y :', LBOUND(V_face_y,1), UBOUND(V_face_y,1), LBOUND(V_face_y,2), UBOUND(V_face_y,2)
    PRINT *, 'Dimensions de C :', LBOUND(C,1), UBOUND(C,1), LBOUND(C,2), UBOUND(C,2)
    PRINT *, 'Dimension de x_noeud :', LBOUND(x_noeud,1), UBOUND(x_noeud,1), LBOUND(x_noeud,2), UBOUND(x_noeud,2)
    PRINT *, 'Dimension de y_noeud :', LBOUND(y_noeud,1), UBOUND(y_noeud,1), LBOUND(y_noeud,2), UBOUND(y_noeud,2)

    CALL VTSWriter(t, 0, get%Nx, get%Ny, x_noeud, y_noeud, C, U_face_x, V_face_y, 'ini')

    DO WHILE (t < time%tf .AND. n < get%Nt)

        CALL STEP_TRANSPORT(recup, get, dim, time, C, C_new, x_reg, y_reg)

        C = C_new
        t = t + dt
        n = n + 1

        CALL VTSWriter(t, n, get%Nx, get%Ny, x_noeud, y_noeud, C, U_face_x, V_face_y, 'int')
        ! --- Sauvegarde à chaque temps ---
        ! WRITE(20,*) "TIME", t

        ! DO j = 1, get%Ny
        !    DO i = 1, get%Nx
        !         WRITE(20,*) x_reg(i), y_reg(j), C(i,j)
        !     END DO
        ! END DO

    END DO

    ! CLOSE(20)

    PRINT *, 'dt =', dt
    PRINT *, 'Simulation terminée après', n, 'pas de temps'
    PRINT *, 'Temps final atteint :', t
    PRINT *, 'Max C =', MAXVAL(C)
    PRINT *, 'Min C =', MINVAL(C)
    
    CALL VTSWriter(t, n, get%Nx, get%Ny, x_noeud, y_noeud, C, U_face_x, V_face_y, 'end')
    DEALLOCATE(x_reg, y_reg, x_noeud, y_noeud, U, U_face_x, V_face_y, C, C_new)

    PRINT *, 'Fin du programme.'
END PROGRAM main_total
