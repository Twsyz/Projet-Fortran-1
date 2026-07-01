! ******************
!  rm -r *.vts
! version paraview: 5.11.2
subroutine VTSWriter(Time,Step,nx,ny,x,y,T,U,V,opt)
!----------------------------------------------------------------------------: Chez Nous#
!  Time    : Reel, temps physique                                            : temps à l'itération
!  Step    : Entier, pas de temps = numero dans le nom de fichier,           : itération 
!            si Step<0 on ecrit                                              : #
!            dans sol_exacte.vts, entier                                     : #
!  nx      : Entier, nombre des cellules en direction x                      : Nx#
!  ny      : Entier, nombre des cellules en direction y                      : Ny#
!  x       : Tableau reel (de taille nx+1,ny+1) des abscisses des noeuds     : x_noeud#
!            des  volumes                                                    :#
!  y       : Tableau reel (de taille nx+1,ny+1) des ordonnees des noeuds     : y_noeud#
!            des  volumes                                                    : #
!  T       : Tableaux reel (de taille nx par ny) des valeurs a tracer        : Concentration au milieu
!            (valeurs au centre des volumes de controle)                     : au milieu des cases
!  U       : Tableaux reel (de taille nx+1   par ny) des valeurs a tracer    :  vitesse x
!            (valeurs au centre des facettes de normale + ou -x              : 
!  V       : Tableaux reel (de taille nx par ny+1  ) des valeurs a tracer    :  vitesse y
!            (valeurs au centre des facettes de normale + ou -y              : #
!  opt     : Variable de type chaine des characteres qui doit prendre        : Booléen qui dit
!            l'une des valeurs suivantes :                                   : quand on prend
!              - 'ini' pour le premier appel a VTSWriter                     : les valeurs ou non
!              - 'int' pour un appel standard a VTSWriter                    : #
!              - 'end' pour le dernier appel a VTSWriter                     : #
!-----------------------------------------------------------------------------#
  implicit none

  real, intent(in)                       :: Time
  integer, intent(in)                    :: Step, nx, ny
  real, dimension(nx+1,ny+1), intent(in) :: x, y
  real, dimension(nx+1,ny)  , intent(in) :: U
  real, dimension(nx,ny+1)  , intent(in) :: V
  real, dimension(nx,ny)    , intent(in) :: T
  character(3), intent(in)               :: opt

  character(100) :: num2char
  character(200) :: FileName, formatperso
  integer :: i, j

  !  --- Ecriture d un fichier temporel au format paraview  ---
  write(num2char,'(i9.9)') Step
  FileName = 'sol_'//trim(num2char)//'.vts'
  open(8,file=FileName)
  write(num2char,*) 3*(nx+1)*(ny+1)
  formatperso = '('//trim(num2char)//'(E15.9,1x))'
  write(8,'(a)') '<?xml version="1.0"?>'
  write(8,'(a)') '<VTKFile type="StructuredGrid">'
  write(8,'(a,6i6,a)') '<StructuredGrid WholeExtent="', 0,nx,0,ny,0,0,'">'
  write(8,'(a,6i6,a)') '<Piece Extent="',0,nx,0,ny,0,0,'">'
  write(8,'(a)') '<Points>'
  write(8,'(a)') '<DataArray type="Float32" NumberOfComponents="3"/>'
  ! Ecriture des coordonnees du maillage
  DO j=1,ny+1
     write(8,formatperso) (x(i,j),y(i,j),0.,i=1,nx+1)
  END DO
  write(8,'(a)') '</Points>'
  write(8,'(a)') '<CellData Scalars="Temperature, U, V">'

  ! Ecriture du scalaire (temperature / concentration)
  write(8,'(a)') '<DataArray type="Float32" Name="Temp, K"/>'
  write(num2char,*) (nx)*(ny)
  DO j=1,ny
     write(8,formatperso) (T(i,j),i=1,nx)
  END DO

  ! Ecriture de la composante U de la vitesse
  write(8,'(a)') '<DataArray type="Float32" Name="Vitesse u, m/s"/>'
  write(num2char,*) (nx)*(ny)
  DO j=1,ny
     write(8,formatperso) ((u(i+1,j)+u(i,j))/2,i=1,nx)
  END DO

  ! Ecriture de la composante V de la vitesse
  write(8,'(a)') '<DataArray type="Float32" Name="Vitesse v, m/s"/>'
  write(num2char,*) (nx)*(ny)
  DO j=1,ny
     write(8,formatperso) ((v(i,j+1)+v(i,j))/2,i=1,nx)
  END DO

  write(8,'(a)') '</CellData>'
  write(8,'(a)') '</Piece>'
  write(8,'(a)') '</StructuredGrid>'
  write(8,'(a)') '</VTKFile>'
  close(8)

  ! - Remplissage du fichier "Collection" determinant l evolution temporelle -
  if (opt == 'ini' ) then
    open(10,file='sol.pvd')
    write(10,'(a)') '<?xml version="1.0"?>'
    write(10,*) '<VTKFile type="Collection" version="0.1" format="ascii">'
    write(10,*) '<Collection>'
  else
    open(10,file='sol.pvd',position='append')
  end if
  if (Step >= 0) write(10,*) '<DataSet timestep="',Time,'" group="" part="0" file="',trim(FileName),'"/>'
  if ( opt == 'end') then
    write(10,*) '</Collection>'
    write(10,*) '</VTKFile>'
  end if
  close(10)

  PRINT *, 'Ecriture du fichier : ', FileName

end subroutine VTSWriter
