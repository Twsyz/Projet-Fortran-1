import numpy as np
import matplotlib.pyplot as plt
import os

print(os.getcwd())
fichier = r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\2\Conv_maillage.dat"
TIME, N, C, Y= [], [], [], []

with open(fichier, 'r') as f:
    lignes = f.readlines()

i = 0
while i < len(lignes):
    if 'TIME' in lignes[i]:
        N.append(float(lignes[i-1].split()[-1]))
        TIME.append(float(lignes[i].split()[-1]))
        i += 1
        
        y_temp, C_temp = [], []
        
        while i < len(lignes) and 'TIME' not in lignes[i] and 'Nx' not in lignes[i]:
            ligne = lignes[i].strip()
            if ligne:  
                data = ligne.split()
                if len(data) >= 2:  
                    y_temp.append(float(data[0]))
                    C_temp.append(float(data[1]))
            i += 1
            
        if y_temp:  
            Y.append(np.array(y_temp))
            C.append(np.array(C_temp))
    else:
        i += 1

print(f"{len(TIME)} maillages différents trouvés.")
print("Nombre de points pour chaque maillage:", N)


plt.figure(figsize=(10, 6))
for j in range(len(TIME)):  # Afficher les 3 premiers jeux de données
    print (Y[j])
    plt.plot(Y[j], C[j], label=f'TIME = {TIME[j]:.2e} s')
    plt.xlabel('Position')
    plt.ylabel('Concentration')
    plt.legend()
    plt.grid(True)
plt.title('Évolution de la concentration au centre pour différents Pe (alpha qui varie)')
plt.savefig('Concentration_vs_Temps_alpha.png')
plt.show()

