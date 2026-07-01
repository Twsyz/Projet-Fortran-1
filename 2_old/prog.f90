! cd /c/Users/moi/Desktop/buro/N7/Methodes\ Volumes\ Finis/BE/2
PROGRAM test_total
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

    REAL :: Cinit, MIN_VALUE, diff_concentration

    INTEGER :: i, j, n, p, Nxx, Nyy, Nxxx, Nyyy 
    REAL :: t, dt, error, i_mid, Pe

    ! -----------------------------
    ! Lecture des données
    ! -----------------------------
    CALL read_data(recup, get, dim, time)
    Nxx= get%Nx
    Nyy= get%Ny
    ! -----------------------------
    ! Allocation
    ! -----------------------------

    PRINT *, 'Allocations réussies'
    PRINT *, 'get%Nx =', get%Nx, 'get%Ny =', get%Ny
    

    ! -----------------------------
    ! Début boucle totale (K va de 1.0e-1 à 1.0e-9)
    ! -----------------------------
    OPEN(30, FILE="Conv_maillage.dat", STATUS="REPLACE")


    ! -----------------------------
    ! Boucle convergence maillage pour Nx et Ny
    ! -----------------------------
    DO p=1,5
        error=0.0

        CALL read_data(recup, get, dim, time) 
        Nxx = get%Nx*p
        Nyy = get%Ny*p
        ALLOCATE(x_reg(Nxx))
        ALLOCATE(y_reg(Nyy))
        ALLOCATE(x_noeud(Nxx+1, Nyy+1))
        ALLOCATE(y_noeud(Nxx+1, Nyy+1))
        ALLOCATE(U(Nxx, Nyy, 2))
        ALLOCATE(U_face_x(Nxx+1, Nyy))
        ALLOCATE(V_face_y(Nxx, Nyy+1))
        ALLOCATE(C(Nxx,Nyy))
        ALLOCATE(C_new(Nxx,Nyy))

        Print *, '--------------------------------------'
        Print *, 'Nx =', Nxx, 'Ny =', Nyy
        WRITE(30,*) "Nx", Nxx, "Ny", Nyy

        ! -----------------------------
        ! Construction du maillage
        ! -----------------------------
        CALL mesh(x_reg, y_reg, dim, get, Nxx, Nyy)
        PRINT *, LBOUND(x_reg,1), UBOUND(x_reg,1)
        PRINT *, LBOUND(y_reg,1), UBOUND(y_reg,1)
        ! CALL noeud(x_noeud, y_noeud, dim, get) 
        ! -----------------------------
        ! Initialisation de la concentration
        ! -----------------------------
        DO j = 1, Nyy
            DO i = 1, Nxx
                C(i,j) = Cinit(recup, get, dim, x_reg(i), y_reg(j), Nxx, Nyy)
            END DO
        END DO
    
        ! -----------------------------
        ! Calcul final
        ! -----------------------------
        CALL Time_Step(recup, dim, get, time, x_reg, y_reg, dt, recup%K)
        PRINT *, "dt =", dt

        t = 0.0
        n = 0   
        error = 1.0

        ! CALL Champ_Vitesse_Faces(recup, dim, get, x_noeud, y_noeud, U_face_x, V_face_y)

        ! PRINT *, 'Dimensions de U_face_x :', LBOUND(U_face_x,1), UBOUND(U_face_x,1), LBOUND(U_face_x,2), UBOUND(U_face_x,2)
        ! PRINT *, 'Dimensions de V_face_y :', LBOUND(V_face_y,1), UBOUND(V_face_y,1), LBOUND(V_face_y,2), UBOUND(V_face_y,2)
        ! PRINT *, 'Dimensions de C :', LBOUND(C,1), UBOUND(C,1), LBOUND(C,2), UBOUND(C,2)
        ! PRINT *, 'Dimension de x_noeud :', LBOUND(x_noeud,1), UBOUND(x_noeud,1), LBOUND(x_noeud,2), UBOUND(x_noeud,2)
        ! PRINT *, 'Dimension de y_noeud :', LBOUND(y_noeud,1), UBOUND(y_noeud,1), LBOUND(y_noeud,2), UBOUND(y_noeud,2)

        !CALL VTSWriter(t, 0, get%Nx, get%Ny, x_noeud, y_noeud, C, U_face_x, V_face_y, 'ini')

        DO WHILE (error > 1.0e-6 .OR. n < 10000)

            CALL STEP_TRANSPORT(recup, get, dim, time, C, C_new, x_reg, y_reg)
            CALL Time_Step(recup, dim, get, time, x_reg, y_reg, dt, recup%K)

            IF (MOD(n, 100) == 0) THEN
                PRINT *, 'erreur =', error
            END IF
            error = diff_concentration(get, C_new, C)

            C = C_new
            t = t + dt
            n = n + 1

            ! CALL VTSWriter(t, n, get%Nx, get%Ny, x_noeud, y_noeud, C, U_face_x, V_face_y, 'int')

            IF (t >= 5.0e-6 .AND. t < 5.0e-6 + dt) THEN
                PRINT *, "Sauvegarde de la concentration à t =", t
                ! --- Sauvegarde de la concentration maximale à chaque temps
                i_mid = Nxx/4
                WRITE(30,*) "TIME", t
                DO j = 1, Nyy
                    WRITE(30,*) y_reg(j), C_new(NINT(i_mid),j)
                END DO
            END IF   

        END DO

        PRINT *, 'dt =', dt
        PRINT *, 'Simulation terminée après', n, 'pas de temps'
        PRINT *, 'Temps final atteint :', t
        PRINT *, 'Max C =', MAXVAL(C)
        PRINT *, 'Min C =', MINVAL(C)

        PRINT *, 'Fin du programme pour Nx =', Nxx, 'Ny =', Nyy
        DEALLOCATE(x_reg, y_reg, x_noeud, y_noeud, U, U_face_x, V_face_y, C, C_new)
    END DO ! Fin de la boucle sur Nx,Ny
    CLOSE(30)

    ! CALL VTSWriter(t, n, get%Nx, get%Ny, x_noeud, y_noeud, C, U_face_x, V_face_y, 'end')

    PRINT *, 'Fin totale'
END PROGRAM test_total
