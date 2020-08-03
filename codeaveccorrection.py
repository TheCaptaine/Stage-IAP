import numpy as np
import matplotlib.pyplot as plt

h = 6.62607004e-34
Kb = 1.38064852e-23
Conversion_Parsec_en_m = 3.086e16
D = 40e6 * Conversion_Parsec_en_m
Cstefan = 5.670374e-8
c = 3e8

choix = 3
Tc1 = 2700 # Tc polaire
Tc2 = 1700 # Tc equateur

R_func = lambda l, Tc : np.sqrt(l/(4*np.pi*Cstefan*Tc**4))
Fflux_func = lambda R, nu_, Tc : (R/D)**2*(2*np.pi*h*nu_**3/c**2*1/(np.exp(h*nu_/(Kb*Tc))-1))
Fmagnitude_func = lambda f : -2.5*np.log10(f)-56.1

temps = np.loadtxt("flux_lambda2.2e-06_m1e+29_v4.47e+07_k0.365.res").T[0]

with open("namefiles.res", "r") as files:
	namefiles = []
	for line in files:
		namefiles += [line[:-1]]

lambd = [2.2e-6, 1.6e-6, 1.2e-6, 1e-6, 0.9e-6, 0.8e-6, 0.7e-6, 0.55e-6, 0.45e-6, 0.365e-6, 0.27e-6]
lettrecode = ["K", "H", "J", "y", "z", "i", "r", "V", "B", "u", "W2"]
codecouleur = ["darkred", "red", "tomato", "orange", "gold", "yellow", "green", "cyan", "royalblue", "blue", "navy"]

date = [i for i in range(57982, 58011, 1)]

data = np.loadtxt('dataa.res')
data2 = np.genfromtxt('dataa2.res', dtype='str')
newdata1 = []
newdata2 = []
newdata3 = []
newdata4 = []
for k in range(len(data[2])):
	if data[2][k] != 0 and data2[k] in lettrecode:
		for j in range(len(date)-1):
			lettre = str(data[0][k])
			taille = len(lettre)
			for t in range(taille):
				if lettre[t] == ".":
					datestr = "0" + str(data[0][k])[t-taille:]
			if data[0][k] >= date[j] and data[0][k] < date[j+1]:
				newdata1 += [j+float(datestr)]
		newdata2 += [data[1][k]]
		newdata3 += [data[2][k]]
		newdata4 += [data2[k]]
alll1 = [[], [], [], [], [], [], [], [], [], [], []]
alll2 = [[], [], [], [], [], [], [], [], [], [], []]
alll3 = [[], [], [], [], [], [], [], [], [], [], []]
for j in range(len(alll1)):
	for k in range(len(newdata1)):
		if newdata4[k] == lettrecode[j]:
			alll1[j] += [newdata1[k]]
			alll2[j] += [newdata2[k]]
			alll3[j] += [newdata3[k]]

#Solo
if choix == 1: #correction partie polaire
	figure1 = plt.figure()
	for k in range(11):
		nu = c/lambd[k]
		luminosite = np.loadtxt(namefiles[int(len(namefiles)/2)+k*4]).T[1]
		T = np.loadtxt(namefiles[int(len(namefiles)/2)+k*4+1]).T[1]
		flux1 = np.loadtxt(namefiles[int(len(namefiles)/2)+k*4+2]).T[1]
		i = 0
		for valeurT in T:
			if valeurT <= Tc1:
				flux1[i] = Fflux_func(R_func(luminosite[i], Tc1), nu, Tc1)
			i+=1
		flux2 = np.loadtxt(namefiles[2+k*4]).T[1]
		plt.plot(temps, Fmagnitude_func(flux2+flux1), color = codecouleur[k], label = lettrecode[k])
		plt.scatter(alll1[k], alll2[k], color = codecouleur[k], zorder = 2)
		plt.errorbar(alll1[k], alll2[k], yerr = alll3[k], fmt = 'none', capsize = 10, ecolor = codecouleur[k], zorder = 1)
elif choix == 2: #correction partie equatoriale
	figure1 = plt.figure()
	for k in range(11):
		nu = c/lambd[k]
		luminosite = np.loadtxt(namefiles[k*4]).T[1]
		T = np.loadtxt(namefiles[k*4+1]).T[1]
		flux1 = np.loadtxt(namefiles[2+k*4]).T[1]
		i = 0
		for valeurT in T:
			if valeurT <= Tc2:
				flux1[i] = Fflux_func(R_func(luminosite[i], Tc2), nu, Tc2)
			i+=1
		flux2 = np.loadtxt(namefiles[int(len(namefiles)/2)+k*4+2]).T[1]
		plt.plot(temps, Fmagnitude_func(flux2+flux1), color = codecouleur[k], label = lettrecode[k])
		plt.scatter(alll1[k], alll2[k], color = codecouleur[k], zorder = 2)
		plt.errorbar(alll1[k], alll2[k], yerr = alll3[k], fmt = 'none', capsize = 10, ecolor = codecouleur[k], zorder = 1)

#En equipe
else:
	figure1 = plt.figure()
	for k in range(11):
		nu = c/lambd[k]
		luminosite = np.loadtxt(namefiles[k*4]).T[1]
		T = np.loadtxt(namefiles[k*4+1]).T[1]
		flux1 = np.loadtxt(namefiles[2+k*4]).T[1]
		i = 0
		for valeurT in T:
			if valeurT <= Tc2:
				flux1[i] = Fflux_func(R_func(luminosite[i], Tc2), nu, Tc2)
			i+=1
		luminosite = np.loadtxt(namefiles[int(len(namefiles)/2)+k*4]).T[1]
		T = np.loadtxt(namefiles[int(len(namefiles)/2)+k*4+1]).T[1]
		flux2 = np.loadtxt(namefiles[int(len(namefiles)/2)+k*4+2]).T[1]
		i = 0
		for valeurT in T:
			if valeurT <= Tc1:
				flux2[i] = Fflux_func(R_func(luminosite[i], Tc1), nu, Tc1)
			i+=1
		plt.plot(temps, Fmagnitude_func(flux2+flux1), color = codecouleur[k], label = lettrecode[k])
		plt.scatter(alll1[k], alll2[k], color = codecouleur[k], zorder = 2)
		plt.errorbar(alll1[k], alll2[k], yerr = alll3[k], fmt = 'none', capsize = 10, ecolor = codecouleur[k], zorder = 1)

plt.legend(loc = 'upper right')
plt.axis([0, 30, 28, 16])
#plt.title("Magnitude d'une kilonova pour différente longueur d'onde en fonction du temps", fontsize = 11)
plt.xlabel("temps [jours]")
plt.ylabel("magnitude")
plt.savefig('figure1.pdf')
plt.show()
