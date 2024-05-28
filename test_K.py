

from function_couplé3 import *
import os


def list_A(N,T):
    nom='assets/stock_data/' + 'QTB_' +str(N) + '_' +str(T) +'.txt'
    if os.path.exists(nom):
        A=[]
        with open(nom, 'r') as fichier:
            lignes = fichier.readlines()
            for ligne in lignes:
                A.append(float(ligne.strip()))  # Convertir la ligne en flottant et l'ajouter à la A

        return A
    else: print("error")


def integrate_K(K):
    result=0
    for i in range(len(K)):
        result+=abs(K[i])

    return result*dt

A=list_A(pas,T)
stop=int(delta_T*5/dt)
K=init_K(stop)

print(integrate_K(K))
print(gamma/integrate_K(K))



print(K[0])
print(K)
# plt.plot(K)
# plt.show()
simulation_QTB(1,Vpprime,A,0,K,stop)
plt.plot(atome[0].x1)
plt.plot(atome[0].x2)
plt.show()






