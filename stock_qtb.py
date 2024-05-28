from qtb import *
import os



def stock_A(A,N):
    nom='assets/stock_data/' + 'QTB_' +str(N) + '_' +str(T) +  '.txt'
    with open(nom,'w') as fichier:
        for i in A:
            fichier.write(str(i) + '\n')



def list_A(N,T):
    nom='assets/stock_data/' + 'QTB_' +str(N) +'_' +str(T) +'.txt'
    if os.path.exists(nom):
        A=[]
        with open(nom, 'r') as fichier:
            lignes = fichier.readlines()
            for ligne in lignes:
                A.append(float(ligne.strip()))  # Convertir la ligne en flottant et l'ajouter à la A
        print(len(A))

        return A
    else:
        A=F_qtb(T,"q")
        stock_A(A,N)
        print("done")
        return A



A=list_A(pas,T)
