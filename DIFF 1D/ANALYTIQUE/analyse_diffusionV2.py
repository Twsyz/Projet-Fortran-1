import numpy as np
import matplotlib.pyplot as plt
import os

L = 1e-3

dossier = r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\DIFF 1D\ANALYTIQUE"

Fichier_tot = {
    "Analyse_R=0.01.dat": 'R=0.01',
    "Analyse_R=0.1.dat": 'R=0.1',
    "Analyse_R=0.3.dat": 'R=0.3',
    'Analyse_R=0.5.dat': 'R=0.5',
    'Analyse_R=0.8.dat': 'R=0.8',
    'Analyse_R=1.dat':   'R=1.0',
    'Analyse_R=1.1.dat': 'R=1.1',
    'Analyse_R=1.15.dat':'R=1.15',
    'Analyse_R=1.25.dat':'R=1.25',
    'Analyse_R=1.5.dat': 'R=1.5',
}

#----------------------------------------------------------------
# Paramètres (si besoin pour l'affichage)
L = 1e-3      # 1 mm

# Correction : utiliser des dictionnaires au lieu de tableaux numpy
temps = {}
y_vals = {}
C_num_vals = {}
C_ana_vals = {}

for k in Fichier_tot:
    # Initialiser les listes pour chaque fichier
    temps[k] = []
    y_vals[k] = []
    C_num_vals[k] = []
    C_ana_vals[k] = []
    
    path = os.path.join(dossier, k)
    with open(path, 'r') as f:
        lignes = f.readlines()

    i = 0
    while i < len(lignes):
        if 'TIME' in lignes[i]:
            t = float(lignes[i].split()[-1])
            i += 1
            y_temp, num_temp, ana_temp = [], [], []
            while i < len(lignes) and lignes[i].strip() and 'TIME' not in lignes[i]:
                data = lignes[i].strip().split()
                if len(data) >= 3:
                    y_temp.append(float(data[0]))
                    num_temp.append(float(data[1]))
                    ana_temp.append(float(data[2]))
                i += 1
            if y_temp:
                temps[k].append(t)
                y_vals[k].append(np.array(y_temp))
                C_num_vals[k].append(np.array(num_temp))
                C_ana_vals[k].append(np.array(ana_temp))
        else:
            i += 1

    print(f"{len(temps[k])} pas de temps lus pour {k}")

    if len(temps[k]) > 0:  # Vérifier qu'il y a des données
        indices = np.linspace(0, len(temps[k])-1, min(4, len(temps[k])), dtype=int)
        C_milieu_num = [profil[len(profil)//2] for profil in C_num_vals[k]]
        C_milieu_ana = [profil[len(profil)//2] for profil in C_ana_vals[k]]

        # --- 4. Evolution temporelle au centre et comparaison à différents temps ---
        plt.figure(figsize=(10, 6))
        
        plt.subplot(2,1,1)
        plt.semilogx(temps[k], C_milieu_num, 'bo-', label='Numérique', markersize=5)
        plt.semilogx(temps[k], C_milieu_ana, 'r--', label='Analytique', linewidth=2)
        plt.xlabel('Temps (s)')
        plt.ylabel(f'C à x = {L*1000:.1f} mm')
        plt.title(f'Évolution au centre - {Fichier_tot[k]}')
        plt.grid(True)
        plt.legend()

        plt.subplot(2,1,2)
        for idx in indices:
            y_mm = y_vals[k][idx]*1000
            plt.plot(y_mm, C_num_vals[k][idx], 'bo-', label=f'Num t={temps[k][idx]:.1e}s', markersize=4)
            plt.plot(y_mm, C_ana_vals[k][idx], 'r--', label=f'Ana t={temps[k][idx]:.1e}s', linewidth=2)
        plt.xlabel('x (mm)')
        plt.ylabel('Concentration C')
        plt.title(f'Comparaison profils - {Fichier_tot[k]}')
        plt.grid(True)
        plt.legend()
        
        plt.tight_layout()
        
        # Sauvegarder avec un nom unique pour chaque fichier
        nom_fichier = f'evolution_comparaison_{Fichier_tot[k].replace("=","")}.png'
        plt.savefig(nom_fichier)
        plt.show()
        plt.close()
    else:
        print(f"Attention : pas de données pour {k}")
print("Fin de l'analyse pour tous les fichiers.")