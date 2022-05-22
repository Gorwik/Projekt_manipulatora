import numpy as n
from matplotlib import pyplot as plt
import scipy.optimize

y,x, = [],[]
v,vy,vx = [],[],[]
ni = []
Fs,R2max,R1n_max,Mgmax = 80.13,0,0,0
l1,l2,l4,l3,l5 = 80,70,90,106.5,14
lin = n.linspace(120,130,1000)

#Obliczenie parametrów niezmiennych w czasie układu
z = n.sqrt(l3**2 + l4**2 - 2*l3*l4*n.cos(n.deg2rad(130)))
ksi = n.arccos((l3**2 + z**2 - l4**2)/(2*l3*z))
q = n.sqrt(z**2 + l5**2 - 2*z*l5*n.cos(n.deg2rad(40))+ksi)
omega = n.arccos((z**2 + q**2 - l5**2)/(2*q*z))

for a in lin:
    p = (l1+l2+a)/2
    P = n.sqrt(p*(p-l1)*(p-l2)*(p-a))
    R = l1*l2*a/(4*P)

    gamma = n.arcsin(l1/(2*R))
    gammadeg = n.rad2deg(gamma)

    alpha = n.arcsin(l2/(2*R))
    alphadeg = n.rad2deg(alpha)

    beta = n.pi-alpha-gamma
    betadeg = 180 - alphadeg - gammadeg
    # print(betadeg)

    delta = n.pi/2 - gamma
    deltadeg = 90 - gammadeg

    fi = n.pi/2 - delta
    fideg = 90 - deltadeg

    psi = delta - ksi
    lambd = psi - omega

    y.append(n.sin(lambd)*q)
    x.append(n.cos(lambd)*q)

    R1 = Fs / (2 * n.cos(alpha))
    Fch = (R1*l2*n.sin(beta))/(q*n.cos(lambd))
    ni.append(Fch*100/Fs)

    # poszukiwanie najwiekszej sily reakcji w podporze nieprzesuwnej
    R2 = n.sqrt(R1**2 + 32.05**2 - 2*R1*32.05*n.cos(n.pi+gamma-beta))
    R2max = R2 if R2 > R2max else R2max

    # obliczenia maksymalnego momentu gnacego i siły normalnej jej towarzyszącej
    Mg = R1 * l2 * n.sin(beta) / 1000
    if Mg > Mgmax:
        Mgmax = Mg
        R1n_max = R1*n.cos(beta)

for p in range(len(y)-1):
    delta = y[p+1] - y[p]
    vy.append(delta/(lin[p+1]-lin[p]))
    delta = x[p + 1] - x[p]
    vx.append(delta / (lin[p + 1] - lin[p]))
    v.append(n.sqrt(vx[p] ** 2 + vy[p] ** 2))

# obliczenie minimalnej wymaganej srednicy tłoka
Fch = 18.75
Fs = Fch*100/(min(ni))
dt = 1000*n.sqrt(4*Fs/(6*10**5*n.pi))
print("Minimalna srednica tłoka siłownika:", dt,"[mm]")

# obliczenie sily fch po doborze siłownika o srednicy (16)
Fchmax = 6*10**5*n.pi*0.016**2*max(ni)/4/100
print("Siła Fch po doborze siłownika o średnicy 16 mm.",Fchmax)

# obliczanie minimalnej srednicy sworznia najbardziej obciazonego z warunku na ścinanie
print("R2 max:",R2max)
Re = 185 * 10**6
kt = 0.23*Re
ds = 1000*n.sqrt(4*R2max/(n.pi*kt))
print("Minimalna srednica sworznia:", ds,"[mm]")

# obliczenia maksymalnych naprężeń
Re = 185
sigma = 0.4*Re*10**6
def kolo(d):
    R = sigma*d**3 - (4*R1n_max/n.pi)*d - 32*Mgmax/n.pi
    return R

def kwadrat(a):
    R = sigma*a**3 - R1n_max*a - 6*Mgmax
    return R

dr = float(1000*scipy.optimize.fsolve(kolo,n.array(1)))
print("Minimalna srednica przekroju poprzecznego ramienia:", dr, "[mm]")

a = float(1000*scipy.optimize.fsolve(kwadrat,n.array(1)))
print("Minimalna dlugosc 'a' przekroju poprzecznego ramienia:", a, "[mm]")

i = n.linspace(0,10,1000)
os_x_do_wykresu = n.linspace(-120,-130,1000)


# plt.subplot(2,2,1)
# plt.title("Wykres przemieszczenia od przesunięcia napędu")
# plt.plot(os_x_do_wykresu,x,os_x_do_wykresu,y)
# plt.legend(["x","y"])
# plt.ylabel("Przemieszczenie chwytaka [mm]")
# plt.xlabel("Przemieszczenie siłownika [mm]")
# plt.grid()
# plt.show()

# plt.subplot(2,2,2)
# plt.title("Wykres stosunku prędkości chwytaka do prędkości napędu")
# plt.plot(os_x_do_wykresu[1:],v,color= "green")
# plt.plot(os_x_do_wykresu[1:],vx)
# plt.plot(os_x_do_wykresu[1:],vy)
# plt.legend(['V',"Vx","Vy"])
# plt.xlabel("Przesunięcie siłownika [mm]")
# plt.ylabel("Stosunek prędkości chwytaka do prędkości napędu [-]")
# plt.grid()
# plt.show()

# plt.subplot(2,2,3)
# plt.title("Wykres stosunku siły na chwytaku do siły napędu")
# plt.plot(os_x_do_wykresu,ni)
# plt.xlabel("Przemieszczenie napędu [mm]")
# plt.ylabel("Stosunek siły na chwytaku do siły siłownika [%]")
# plt.grid()
# plt.show()