import matplotlib.pyplot as plt

fitness_same = [2.13E-08, 311.387761409, 0.3747056007 ] 
fitness_modi = [0, 1.62E-10, 0.0429593738]
fitness_maxl = [0, 2.63E-11, 6.18E-12]
fintess_nsteps = [0.064549041, 3756.55380419, 1.68E-10]
time_same = [4.1913487911, 0.0644550323, 0.412115097]
time_modi=[0, 0.282102108, 0.3806059361]
time_maxl=[ 0, 0.5261609554, 0.5216841698]
time_nsteps = [0.4833509922, 0.0189540386, 0.0940279961]
labels_clust = ["Same prameters","Modified parameters", "MaxLocIter 100", "NSteps=PopSize"]
labels_in = ["RMGS_all", "RMGS", "GGS"]

fig, ax = plt.subplots()
ax.semilogy(time_same, fitness_same, 'o')
ax.semilogy(time_modi, fitness_modi, 'ro')
ax.semilogy(time_maxl, fitness_maxl, 'yo')
ax.semilogy(time_nsteps, fintess_nsteps, 'go')
plt.legend(labels_clust)
for i, txt in enumerate(labels_in):
    ax.annotate(txt, (time_same[i],fitness_same[i]))
    ax.annotate(txt, (time_modi[i],fitness_modi[i]))
    ax.annotate(txt, (time_maxl[i],fitness_maxl[i]))
    ax.annotate(txt, (time_nsteps[i],fintess_nsteps[i]))
plt.draw()
plt.show()