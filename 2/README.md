--------------------------------------------------------------------------
# 1.Informations générales
Nom du programme : prog.exe
Date de publication: 06/05/2026
Auteur : JACOB Yvan et AUDEBERT–CASTELLI Mathéo 
Objectif : Programme (en Fortran 90) de résolution d'une équation
d'advection-diffusion 2D. Methode des volumes finis. Schema Euler
en temps, amont pour l'advection, centre pour la diffusion.
Fichier d'entree: "donnees.dat"
Fichiers resultats : "sol_000000.vts", "sol_000001.vts", etc., sol.pvd;
Pour compiler le programme :
Dans un terminal, se placer dans le dossier contenant les sources
du programme (*.f90, Makefile), puis taper:
make clean
make
Pour executer le programme :
Dans un terminal, se placer dans le dossier contenant le fichier
executable (prog.exe) et le fchier d'entree (donnees.dat), puis
taper:
./prog.exe
Pour visualiser les resultats:
Dans un terminal, se placer dans le dossier contenant les
fichiers de resultat (*.vts), puis taper:
paraview &

Dans les fichiers Python, il ne faut pas hésiter à changer les chemins d'accès aux fichiers de données et de résultats pour les adapter à votre environnement de travail.

Pour avoir les animations dans le pdf, veuillez utiliser Adobe Acrobat Reader, car d'autres lecteurs de pdf peuvent ne pas supporter les animations intégrées.
--------------------------------------------------------------------------

--------------------------------------------------------------------------
# 2.Architecture
==========================================================================
## 1. Solution pour l'advéction pure
Dossier DIFF 1D ->
### 1. Selon X -> Dossier X

### 1. Selon Y -> Dossier Y

### 1. Comparaison avec python -> Dossier Analytique
==========================================================================

==========================================================================
## 2. Solution pour la diffusion pure

### 1. Selon X -> Dossier X

### 1. Selon Y -> Dossier Y
==========================================================================

==========================================================================
## 3. Calculer le champ stationnaire de concentration pour un nombre de Péclet massique donné -> Dossier 1
Calcul du champ stationnaire en vérifiant qu'il n'y a pas de changement de concentration à une erreur ('error') près.
==========================================================================

==========================================================================
## 4. Etude de la convergence en maillage du simulateur -> Dossier 2
Calcul pour plusieurs maillages:
Changer Nx et Ny dans "donnees.dat", renommer le fichier "profil.dat" en "profil20.dat" pour un maillage de 20 par 20. Après ajouter à "analyse.py": 
yX, CX = np.loadtxt(r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\2\profilX.dat", unpack=True)
et rajouter aussi le plot. 
On observe une convergence pour un maillage autour de 100
==========================================================================

==========================================================================
## 5. Campagne de simulations en faisant varier le nombre de Péclet (par variation du coefficient de Kappa) -> Dossier 3
Nous l'avons fait pour plus de nombre de Péclet afin d'observer des résultats pour la comparaison des temps adimensionnels question 5.
Les modifications sont principalement dans le programme principal sont l'ajout d'une boucle pour faire varier le nombre de Péclet et le fichier de sortie "It_PE.dat" pour stocker les résultats.
Utiliser le fichier "analyse_Pe.dat" pour faire les graphiques de la question 3.
==========================================================================

==========================================================================
## 5. Campagne de simulations en faisant varier le nombre de Péclet (par variation du coefficient de vitesse alpha) -> Dossier 4
Nous l'avons fait pour plus de nombre de Péclet afin d'observer des résultats pour la comparaison des temps adimensionnels question 5. Le fichie de sortie est "It_PE_A.dat" ou "It_PE.dat".
Utiliser le fichier "analyse_Pe.dat" pour faire les graphiques de la question 4.
==========================================================================

==========================================================================
## 5. Campagne de simulations en faisant varier le nombre de Péclet (par variation du coefficient de vitesse alpha) -> Dossier 5
Exhibiton des t95%, comparaison des temps adimensionnels pour les différentes méthodes de calcul du nombre de Péclet.
==========================================================================
--------------------------------------------------------------------------