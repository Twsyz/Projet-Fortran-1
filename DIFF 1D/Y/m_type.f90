MODULE m_type
	IMPLICIT NONE
	
	TYPE taille
		REAL :: L, D
	END TYPE taille

	TYPE phys
		REAL ::  C0, C1, K, alpha
		REAL, DIMENSION(2) :: U
	END TYPE phys

	TYPE phys_time
		REAL :: tf, CFL, R
	END TYPE phys_time

	TYPE num
		INTEGER :: Nx,Ny,Nt
	END TYPE num

	TYPE concentration
		REAL :: CA, CB, CC, CD, CE
	END TYPE concentration

	TYPE flux_diff
		REAL :: Q1, Q2, Q3, Q4
	END TYPE flux_diff

	TYPE flux_adv
		REAL :: F1, F2, F3, F4
	END TYPE flux_adv

END MODULE m_type









