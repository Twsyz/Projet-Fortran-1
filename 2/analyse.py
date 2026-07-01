import numpy as np
import matplotlib.pyplot as plt
import os

print(os.getcwd())

y1, C1 = np.loadtxt(r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\2_save\profil10.dat", unpack=True)
y2, C2 = np.loadtxt(r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\2_save\profil20.dat", unpack=True)
y3, C3 = np.loadtxt(r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\2_save\profil40.dat", unpack=True)
y4, C4 = np.loadtxt(r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\2_save\profil80.dat", unpack=True)
y5, C5 = np.loadtxt(r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\2_save\profil100.dat", unpack=True)
y6, C6 = np.loadtxt(r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\2_save\profil120.dat", unpack=True)
y7, C7 = np.loadtxt(r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\2_save\profil200.dat", unpack=True)
y8, C8 = np.loadtxt(r"c:\Users\moi\Desktop\buro\N7\Methodes Volumes Finis\BE\2_save\profil300.dat", unpack=True)

plt.plot(y1,C1,label="10x10")
plt.plot(y2,C2,label="20x20")
plt.plot(y3,C3,label="40x40")
plt.plot(y4,C4,label="80x80")
plt.plot(y5,C5,label="100x100")
plt.plot(y6,C6,label="120x120")


plt.xlabel("y")
plt.ylabel("C en x=L/2")
plt.title("Concentration en x=L/2, pour l'instant final (convergence en maillage)")
plt.legend()
plt.show()
# plt.savefig("convergence_maillage_1.png")

plt.plot(y4,C4,label="80x80")
plt.plot(y5,C5,label="100x100")
plt.plot(y6,C6,label="120x120")
plt.plot(y7,C7,label="200x200")
plt.plot(y8,C8,label="300x300")


plt.xlabel("y")
plt.ylabel("C en x=L/2")
plt.title("Concentration en x=L/2, pour l'instant final (convergence en maillage)")
plt.legend()
plt.show()
# plt.savefig("convergence_maillage_2.png")


