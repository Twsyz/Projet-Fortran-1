import numpy as np
import matplotlib.pyplot as plt
import os

print(os.getcwd())
fichier = r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\3\It_Pe.dat"
PE, t, C = [], [], []

with open(fichier, 'r') as f:
    lignes = f.readlines()

i = 0
while i < len(lignes):
    if 'Pe' in lignes[i]:
        PE.append(float(lignes[i].split()[-1]))
        i += 1
        
        t_temp, C_temp = [], []
        
        while i < len(lignes) and 'Pe' not in lignes[i]:
            ligne = lignes[i].strip()
            if ligne:  
                data = ligne.split()
                if len(data) >= 2:  
                    t_temp.append(float(data[0]))
                    C_temp.append(float(data[1]))
            i += 1
            
        if t_temp:  
            t.append(np.array(t_temp))
            C.append(np.array(C_temp))
    else:
        i += 1

print(f"{len(t)} Peclet différents")
print("PE =", PE)

# Calcul de t95%
# Recuperer la derniere valeur de C pour chaque Pe
C_final = [C[j][-1] for j in range(len(C))]
# Parcourir les concentrations pour trouver le temps où C atteint 95% de C_final
t95 = []
j=0
for j in range(len(C)):
    C_95 = 0.95 * C_final[j]
    t_95 = None
    for k in range(len(C[j])):
        if C[j][k] >= C_95:
            t_95 = t[j][k]
            break
    t95.append(t_95)
print("t95% pour chaque Pe :", t95)

plt.figure(figsize=(10, 6))
for j in range(len(t)):  # Afficher les 3 premiers jeux de données
    plt.plot(t[j], C[j], label=f'Pe = {PE[j]:.2e}, t95={t95[j]:.2e}s')
    plt.xlabel('Temps')
    plt.ylabel('Concentration')
    plt.legend()
    plt.grid(True)
plt.title('Évolution de la concentration au centre pour différents Pe (Kappa qui varie)')
plt.savefig('Concentration_vs_Temps_Kappa.png')
plt.show()

