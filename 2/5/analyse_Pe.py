import numpy as np
import matplotlib.pyplot as plt
import os

print(os.getcwd())

L=1.0*10**(-3)  # Longueur du domaine
K=1.0*10**(-5)  # Diffusivité
A=1.0           # Vitesse d'advection

fichier_alpha = r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\5\It_Pe_A_V2.dat"
fichier_kappa = r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\5\It_Pe_K_V2.dat"
PE_alpha, t_alpha, C_alpha = [], [], []

with open(fichier_alpha, 'r') as f:
    lignes_alpha = f.readlines()

i = 0
while i < len(lignes_alpha):
    if 'Pe' in lignes_alpha[i]:
        PE_alpha.append(float(lignes_alpha[i].split()[-1]))
        i += 1
        
        t_temp, C_temp = [], []
        
        while i < len(lignes_alpha) and 'Pe' not in lignes_alpha[i]:
            ligne_alpha = lignes_alpha[i].strip()
            if ligne_alpha:  
                data_alpha = ligne_alpha.split()
                if len(data_alpha) >= 2:  
                    t_temp.append(float(data_alpha[0]))
                    C_temp.append(float(data_alpha[1]))
            i += 1
            
        if t_temp:  
            t_alpha.append(np.array(t_temp))
            C_alpha.append(np.array(C_temp))
    else:
        i += 1

PE_kappa, t_kappa, C_kappa = [], [], []

with open(fichier_kappa, 'r') as f:
    lignes_kappa = f.readlines()

i = 0
while i < len(lignes_kappa):
    if 'Pe' in lignes_kappa[i]:
        PE_kappa.append(float(lignes_kappa[i].split()[-1]))
        i += 1
        
        t_temp, C_temp = [], []
        
        while i < len(lignes_kappa) and 'Pe' not in lignes_kappa[i]:
            ligne_kappa = lignes_kappa[i].strip()
            if ligne_kappa:  
                data_kappa = ligne_kappa.split()
                if len(data_kappa) >= 2:  
                    t_temp.append(float(data_kappa[0]))
                    C_temp.append(float(data_kappa[1]))
            i += 1
            
        if t_temp:  
            t_kappa.append(np.array(t_temp))
            C_kappa.append(np.array(C_temp))
    else:
        i += 1

print(f"{len(t_alpha)} Peclet différents (alpha)")
print("PE_alpha =", PE_alpha)
print(f"{len(t_kappa)} Peclet différents (kappa)")
print("PE_kappa =", PE_kappa)


# Calcul de t95%
# Recuperer la derniere valeur de C pour chaque Pe
C_final_kappa = [C_kappa[j][-1] for j in range(len(C_kappa))]
C_final_alpha = [C_alpha[j][-1] for j in range(len(C_alpha))]
# Parcourir les concentrations pour trouver le temps où C atteint 95% de C_final
t95_kappa, t95_alpha = [], []

for j in range(len(C_kappa)):
    C_95 = 0.95 * C_final_kappa[j]
    t_95 = None
    for k in range(len(C_kappa[j])):
        if C_kappa[j][k] >= C_95:
            t_95 = t_kappa[j][k]
            break
    t95_kappa.append(t_95)
print("t95% pour chaque Pe (kappa) :", t95_kappa)

for j in range(len(C_alpha)):
    C_95 = 0.95 * C_final_alpha[j]
    t_95 = None
    for k in range(len(C_alpha[j])):
        if C_alpha[j][k] >= C_95:
            t_95 = t_alpha[j][k]
            break
    t95_alpha.append(t_95)
print("t95% pour chaque Pe (alpha) :", t95_alpha)

t_95_kappa_adim = [t95_kappa[j] * K / L**2 for j in range(len(t95_kappa))]
t_95_alpha_adim = [t95_alpha[j] * A / L for j in range(len(t95_alpha))]

print("t95% adimensionnel pour chaque Pe (diffusion) :", t_95_kappa_adim)
print("t95% adimensionnel pour chaque Pe (advéction) :", t_95_alpha_adim)

plt.figure(figsize=(10, 6))
plt.minorticks_on()
plt.grid(True, which='both', linestyle='--', linewidth=0.5)

plt.loglog(PE_kappa, t_95_kappa_adim, marker='o', label='t95% adimensionnel (kappa)')   
plt.loglog(PE_alpha, t_95_alpha_adim, marker='o', label='t95% adimensionnel (alpha)')
plt.xlabel("Nombre de Peclet (Pe)")
plt.ylabel("t95% adimensionnel")
plt.title("Temps adimensionnels en fonction du nombre de Péclet")
plt.legend()
plt.savefig("temps_adimensionnels_superposition")
plt.show()
