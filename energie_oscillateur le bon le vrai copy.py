


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

print(len(K))
xhi = np.sqrt(hb/(m*w1))

def E(c3,A):
    C3=c3*hb*w1/(xhi**3)
    simulation_QTB(1,Vpprime,A,C3,K,stop)
    print("ok")
    E1 = []
    E2 = []
    E3 = []
    e1=0
    e2=0
    for i in range (len(atome[0].v1)) :
        e1 += 1/(2*xhi**2 ) * (((atome[0].v1[i])**2)/(w1**2) + (atome[0].x1[i])**2)
        e2 += 1/(2*xhi**2 ) * (((atome[0].v2[i])**2)/(w1**2) + W**2*(atome[0].x2[i])**2)  
        E1.append(e1/(i+1))
        E2.append(e2/(i+1))     
        E3.append((e1+e2)/(i+1))      
    print("oui")            
    return [E1,E2,E3]  



p=int(pas/1000)

# A=F_qtb(T,"q")

c=np.linspace(0,140,1)*10**(-5)
print(len(c),c)

L1=[]
L2=[]
L3=[]


for j in range(len(c)):
    print(j,len(c))
    i=c[j]
    result=E(i,A)
    L1.append(result[0][-1])
    L2.append(result[1][-1])
    L3.append(result[2][-1])
    print(c,L1,L2,L3)



#Tracer la fonction

plt.plot(c, L1,label='E1')
plt.plot(c, L2,label='E2')
plt.plot(c, L3,label='E3')
plt.legend()
plt.title("Énergie moeynne des oscillateurs")
plt.xlabel("x")
plt.ylabel("y")
plt.grid(True)
plt.show()
