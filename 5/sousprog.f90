! ============================================================
! Sous-programmes pour calcul total
! ============================================================

! -------------------------------------------------------------------------------
! Subroutine pour lire les données à partir d'un fichier et les stocker dans des types définis
! -------------------------------------------------------------------------------

SUBROUTINE read_data(recup, get, dim, time)
	USE m_type
	IMPLICIT NONE

	TYPE(phys), INTENT(OUT) :: recup
	TYPE(num), INTENT(OUT) :: get
    TYPE(taille), INTENT(OUT) :: dim
	TYPE(phys_time), INTENT(OUT) :: time


	OPEN(10,FILE="donnees.dat")

	READ(10,*) dim%L
	READ(10,*) dim%D
    READ(10,*) recup%C0
    READ(10,*) recup%C1
	READ(10,*) recup%K
	READ(10,*) recup%alpha
	READ(10,*) time%tf
    READ(10,*) time%CFL
    READ(10,*) time%R
	READ(10,*) 
	READ(10,*) get%Nx
	READ(10,*) get%Ny
    READ(10,*) get%Nt

	CLOSE(10)

END SUBROUTINE read_data

! -------------------------------------------------------------------------------
! Subroutine pour générer les maillages en x et y à partir des paramètres de taille et du nombre de points
! -------------------------------------------------------------------------------

SUBROUTINE mesh(x_reg, y_reg, dim, get)
    USE m_type
    IMPLICIT NONE
    
    TYPE(taille), INTENT(IN) :: dim
    TYPE(num), INTENT(IN) :: get
    REAL, DIMENSION(get%Nx), INTENT(OUT) :: x_reg
    REAL, DIMENSION(get%Ny), INTENT(OUT) :: y_reg
    INTEGER :: i,j
    REAL :: dx, dy

    dx = 2.0*dim%L / REAL(get%Nx)
    dy = dim%L / REAL(get%Ny)

    DO i=1,get%Nx
        x_reg(i) = dx/2.0 + (i-1)*dx  
    END DO
    DO j=1,get%Ny
        y_reg(j) = dy/2.0 + (j-1)*dy  
    END DO

END SUBROUTINE mesh

SUBROUTINE noeud(x_noeud, y_noeud, dim, get)
    USE m_type
    IMPLICIT NONE
    
    TYPE(taille), INTENT(IN) :: dim
    TYPE(num), INTENT(IN) :: get

    REAL, DIMENSION(get%Nx+1, get%Ny+1), INTENT(OUT) :: x_noeud, y_noeud
    INTEGER :: i,j

    DO i=1,get%Nx+1
        DO j=1,get%Ny+1
            x_noeud(i,j) = 2.0*dim%L * REAL(i-1)/REAL(get%Nx)
            y_noeud(i,j) = dim%L * REAL(j-1)/REAL(get%Ny)
        END DO
    END DO


END SUBROUTINE

! -------------------------------------------------------------------------------
! Fonction pour calculer la condition initiale en fonction de la position et des paramètres physiques
! -------------------------------------------------------------------------------

FUNCTION Cinit(recup, get, dim, x, y)
	USE m_type 
	IMPLICIT NONE 

	TYPE(phys), INTENT(IN) :: recup
	TYPE(num), INTENT(IN) :: get
	TYPE(taille), INTENT(IN) :: dim

    REAL :: x, y, dist

    REAL :: Cinit
    REAL :: DELTA
    
    DELTA = (2./3.)*min((2*dim%L)/real(get%Nx), (dim%L)/real(get%Ny))

    dist = sqrt((x - dim%L)**2 + (y - dim%L/2.0)**2)

    Cinit = recup%C1 + 0.5*(recup%C0 - recup%C1) * (1.0 + erf((dist - dim%D/2.0) / DELTA))
END FUNCTION Cinit

! -------------------------------------------------------------------------------
! Fonction pour calculer la vitesse à partir des paramètres physiques et de la position
! -------------------------------------------------------------------------------

SUBROUTINE Vitesse(recup, dim, x, y, V)
    USE m_type
    IMPLICIT NONE

    TYPE(phys), INTENT(IN) :: recup
    TYPE(taille), INTENT(IN) :: dim

    REAL, INTENT(IN) :: x, y
    REAL, INTENT(OUT) :: V(2)

    REAL :: pi
    pi = ACOS(-1.0)

    V(1) =  recup%alpha * SIN(pi*x/dim%L) * COS(pi*y/dim%L)
    V(2) = -recup%alpha * COS(pi*x/dim%L) * SIN(pi*y/dim%L)
END SUBROUTINE Vitesse

! -------------------------------------------------------------------------------
! Fonction pour trouver la valeur minimale dans une matrice
! -------------------------------------------------------------------------------
! Fonction plus utilisée...
FUNCTION MIN_VALUE(C, taille_1, taille_2)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: taille_1, taille_2

    REAL, DIMENSION(taille_1, taille_2) :: C
    INTEGER :: m,p
    REAL :: MIN_VALUE
    
    MIN_VALUE = C(1,1)
    DO m=1, taille_1
        DO p=1, taille_2
            IF (C(m,p) < MIN_VALUE) THEN
                MIN_VALUE = C(m,p)
            END IF
        END DO
    END DO
END FUNCTION MIN_VALUE

! -------------------------------------------------------------------------------
! Subroutine pour calculer le champ de vitesse à chaque point du maillage
! -------------------------------------------------------------------------------

SUBROUTINE Champ_Vitesse(recup, dim, get, x_reg, y_reg, U)
    USE m_type
    IMPLICIT NONE
    
    TYPE(phys), INTENT(IN) :: recup
    TYPE(taille), INTENT(IN) :: dim
    TYPE(num), INTENT(IN) :: get

    REAL, DIMENSION(get%Nx), INTENT(IN) :: x_reg
    REAL, DIMENSION(get%Ny), INTENT(IN) :: y_reg
    REAL, DIMENSION(get%Nx,get%Ny,2), INTENT(OUT) :: U

    INTEGER :: i,j
    REAL :: V(2)

    DO i=1,get%Nx
        DO j=1,get%Ny
            CALL Vitesse(recup, dim, x_reg(i), y_reg(j), V)
            U(i,j,:) = V
        END DO
    END DO

END SUBROUTINE Champ_Vitesse


! -------------------------------------------------------------------------------
! Subroutine pour calculer les vitesses aux faces du maillage
! (aux faces verticales pour U, aux faces horizontales pour V)
! -------------------------------------------------------------------------------

SUBROUTINE Champ_Vitesse_Faces(recup, dim, get, x_noeud, y_noeud, U_face, V_face)
    USE m_type
    IMPLICIT NONE
    
    TYPE(phys), INTENT(IN) :: recup
    TYPE(taille), INTENT(IN) :: dim
    TYPE(num), INTENT(IN) :: get

    REAL, DIMENSION(get%Nx+1, get%Ny+1), INTENT(IN) :: x_noeud, y_noeud
    REAL, DIMENSION(get%Nx+1, get%Ny), INTENT(OUT) :: U_face  ! Vitesse u aux faces verticales
    REAL, DIMENSION(get%Nx, get%Ny+1), INTENT(OUT) :: V_face  ! Vitesse v aux faces horizontales

    INTEGER :: i, j
    REAL :: x_face, y_face
    REAL :: V(2)



    DO j = 1, get%Ny
        DO i = 1, get%Nx+1
            x_face = x_noeud(i, j)
            y_face = y_noeud(i, j)
            CALL Vitesse(recup, dim, x_face, y_face, V)
            U_face(i, j) = V(1)  ! Composante u (vitesse en x)
        END DO
    END DO

    DO j = 1, get%Ny+1
        DO i = 1, get%Nx

            x_face = x_noeud(i, j)
            y_face = y_noeud(i, j)
            CALL Vitesse(recup, dim, x_face, y_face, V)
            V_face(i, j) = V(2)  ! Composante v (vitesse en y)
        END DO
    END DO

END SUBROUTINE Champ_Vitesse_Faces

! -------------------------------------------------------------------------------
! Calcul du pas de temps en fonction de la vitesse et des paramètres physiques
! -------------------------------------------------------------------------------

! Problème avec CFL et R qui sont dans time et pas dans recup... à revoir
SUBROUTINE Time_Step(recup, dim, get, time, x_reg, y_reg, dt)
    USE m_type
    IMPLICIT NONE
    
    TYPE(phys), INTENT(IN) :: recup
    TYPE(taille), INTENT(IN) :: dim
    TYPE(num), INTENT(IN) :: get
    TYPE(phys_time), INTENT(IN) :: time

    REAL, DIMENSION(get%Nx), INTENT(IN) :: x_reg
    REAL, DIMENSION(get%Ny), INTENT(IN) :: y_reg
    REAL, INTENT(OUT) :: dt

    REAL :: dx, dy
    REAL :: V(2)
    REAL :: dt_local, dt_min
    INTEGER :: i,j

    dx = 2.0*dim%L / REAL(get%Nx-1)
    dy = dim%L / REAL(get%Ny-1)

    dt_min = 1.0E30

    DO i=1,get%Nx
        DO j=1,get%Ny
            CALL Vitesse(recup, dim, x_reg(i), y_reg(j), V)

            dt_local = 1.0 / ( ABS(V(1))/(dx*time%CFL) + &
                               ABS(V(2))/(dy*time%CFL) + &
                               recup%K/((dx*dx)*time%R) + &
                               recup%K/((dy*dy)*time%R) )

               !dt_local = time%CFL / ( ABS(V(1))/dx + ABS(V(2))/dy + &
               !         2.0*recup%K/(dx*dx) + 2.0*recup%K/(dy*dy) )


            IF (dt_local < dt_min) dt_min = dt_local
        END DO
    END DO

    dt = dt_min
END SUBROUTINE Time_Step

! -------------------------------------------------------------------------------
! Subroutine pour  calculer le flux de diffusion en chaque point
! -------------------------------------------------------------------------------

SUBROUTINE DIFFUSION(recup, get, dim, conc, flux_d)
    USE m_type 
    IMPLICIT NONE 
    TYPE(phys), INTENT(IN) :: recup
    TYPE(num), INTENT(IN) :: get
    TYPE(taille), INTENT(IN) :: dim
    TYPE(concentration),INTENT(IN) :: conc
    TYPE(flux_diff), INTENT(OUT) :: flux_d
    
    REAL :: dx,dy

    dx = 2*dim%L / real(get%Nx-1)
    dy = dim%L / real(get%Ny-1)
    ! Attention a qui appartient a qui dans les indices de flux_d et concentration
    flux_d%Q1 = -recup%K * (conc%CE - conc%CA) / dx ! Flux en j de (i-1) -> i OUEST
    flux_d%Q2 = -recup%K * (conc%CB - conc%CE) / dx ! Flux en j de i -> (i+1) EST
    flux_d%Q3 = -recup%K * (conc%CE - conc%CC) / dy ! Flux en i de (j-1) -> j SUD
    flux_d%Q4 = -recup%K * (conc%CD - conc%CE) / dy ! Flux en i de j -> (j+1) NORD

END SUBROUTINE DIFFUSION

! -------------------------------------------------------------------------------
! Fonction pour calculer l'advetion en chaque point du maillage en fonction de la vitesse et de la concentration
! -------------------------------------------------------------------------------

SUBROUTINE ADVECTION(recup, dim, conc, x, y, flux_a)
    USE m_type
    IMPLICIT NONE

    TYPE(phys), INTENT(IN) :: recup
    TYPE(taille), INTENT(IN) :: dim
    TYPE(concentration), INTENT(IN) :: conc
    TYPE(flux_adv), INTENT(OUT) :: flux_a

    REAL, INTENT(IN) :: x, y
    REAL :: V(2)

    ! Récupération vitesse locale
    CALL Vitesse(recup, dim, x, y, V)

    ! =====================
    ! Direction X (Marche pas)
    ! =====================

    !IF (V(1) >= 0.0) THEN
    !    flux_a%F1 =  -V(1) * conc%CA   ! OUEST
    !    flux_a%F2 =  V(1) * conc%CE   ! EST
    !ELSE
    !    flux_a%F1 =  -V(1) * conc%CE
    !    flux_a%F2 =  V(1) * conc%CB
    !END IF

    ! =====================
    ! Direction Y (Marche pas)
    ! =====================

    !IF (V(2) >= 0.0) THEN
    !    flux_a%F3 =  -V(2) * conc%CC   ! SUD
    !    flux_a%F4 =  V(2) * conc%CE   ! NORD
    !ELSE
    !    flux_a%F3 =  -V(2) * conc%CE
    !   flux_a%F4 =  V(2) * conc%CD
    !END IF

    ! =====================
    ! Direction X 
    ! =====================

    IF (V(1) >= 0.0) THEN
        flux_a%F1 =  V(1) * conc%CA   ! OUEST
        flux_a%F2 =  V(1) * conc%CE   ! EST
    ELSE
        flux_a%F1 =  V(1) * conc%CE
        flux_a%F2 =  V(1) * conc%CB
    END IF

    ! =====================
    ! Direction Y (Marche pas)
    ! =====================

    IF (V(2) >= 0.0) THEN
        flux_a%F3 =  V(2) * conc%CC   ! SUD
        flux_a%F4 =  V(2) * conc%CE   ! NORD
    ELSE
        flux_a%F3 =  V(2) * conc%CE
        flux_a%F4 =  V(2) * conc%CD
    END IF

END SUBROUTINE ADVECTION

! -------------------------------------------------------------------------------
! Final
! -------------------------------------------------------------------------------

! Stabilité

! Advection : dt <= CFL * min(dx/|Vx|, dy/|Vy|)
! Diffusion : dt <= 0.5 * min(dx^2, dy^2)/K

SUBROUTINE STEP_TRANSPORT(recup, get, dim, time, C, C_new, x_reg, y_reg)

    USE m_type
    IMPLICIT NONE

    TYPE(phys), INTENT(IN) :: recup
    TYPE(num), INTENT(IN) :: get
    TYPE(taille), INTENT(IN) :: dim
    TYPE(phys_time), INTENT(IN) :: time

    REAL, INTENT(IN) :: x_reg(get%Nx), y_reg(get%Ny)
    REAL, INTENT(IN) :: C(get%Nx, get%Ny)
    REAL, INTENT(OUT) :: C_new(get%Nx, get%Ny)

    TYPE(concentration) :: conc
    TYPE(flux_diff) :: flux_d
    TYPE(flux_adv) :: flux_a

    REAL :: dx, dy, dt
    INTEGER :: i, j

    dx = 2*dim%L / REAL(get%Nx-1)
    dy = dim%L / REAL(get%Ny-1)

    CALL Time_Step(recup, dim, get, time, x_reg, y_reg, dt)

    C_new=C

    DO j = 2, get%Ny-1
        DO i = 2, get%Nx-1

            ! Remplir structure conc
            conc%CE = C(i,j)
            conc%CA = C(i-1,j)
            conc%CB = C(i+1,j)
            conc%CC = C(i,j-1)
            conc%CD = C(i,j+1)

            ! Diffusion
            CALL DIFFUSION(recup, get, dim, conc, flux_d)

            ! Advection
            CALL ADVECTION(recup, dim, conc, x_reg(i), y_reg(j), flux_a)

            ! Mise à jour explicitee
            C_new(i,j) = C(i,j)  &
              + dt/dx * (flux_a%F1 - flux_a%F2) &
              + dt/dy * (flux_a%F3 - flux_a%F4) &
              + dt/dx * (flux_d%Q1 - flux_d%Q2) &
              + dt/dy * (flux_d%Q3 - flux_d%Q4)

            IF (ABS(C_new(i,j)) > 1.0E10) THEN
                    PRINT *, 'Explosion à i=', i, 'j=', j, 'C=', C_new(i,j)
                    PRINT *, 'F1,F2,F3,F4=', flux_a%F1, flux_a%F2, flux_a%F3, flux_a%F4
                    PRINT *, 'Q1,Q2,Q3,Q4=', flux_d%Q1, flux_d%Q2, flux_d%Q3, flux_d%Q4
                    STOP
                END IF
        END DO
    END DO

    ! Bords gauche et droit
    DO j = 2, get%Ny-1
        C_new(1,j) = C_new(2,j)      ! Flux nul en x=0
        C_new(get%Nx,j) = C_new(get%Nx-1,j)  ! Flux nul en x=2L
    END DO
    
    ! Bords haut et bas
    DO i = 2, get%Nx-1
        C_new(i,1) = C_new(i,2)      ! Flux nul en y=0
        C_new(i,get%Ny) = C_new(i,get%Ny-1)  ! Flux nul en y=L
    END DO

    ! Coins
    C_new(1,1) = 0.5*(C_new(2,1) + C_new(1,2))
    C_new(1,get%Ny) = 0.5*(C_new(2,get%Ny) + C_new(1,get%Ny-1))
    C_new(get%Nx,1) = 0.5*(C_new(get%Nx-1,1) + C_new(get%Nx,2))
    C_new(get%Nx,get%Ny) = 0.5*(C_new(get%Nx-1,get%Ny) + C_new(get%Nx,get%Ny-1))


END SUBROUTINE STEP_TRANSPORT

FUNCTION diff_concentration(get, Cnew, Cold)

    USE m_type
    IMPLICIT NONE

    TYPE(num), INTENT(IN) :: get
    REAL, INTENT(IN) :: Cnew(get%Nx,get%Ny), Cold(get%Nx,get%Ny)

    REAL :: diff_concentration, diffmax, diff
    INTEGER :: i, j

diffmax = 1.0e-30

DO j = 1, get%Ny
    DO i = 1, get%Nx

        diff = ABS(Cnew(i,j) - Cold(i,j))
        IF (diff > diffmax) THEN
            diffmax = diff
        END IF

    END DO
END DO
diff_concentration = diffmax
END FUNCTION diff_concentration



