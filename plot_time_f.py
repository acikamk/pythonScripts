import matplotlib.pyplot as plt

fitness_same = [1.2684476800372846, 29.111040935987955, 6.9743849149343404]
time_same = [0.07431412935256958, 0.008248186111450196, 0.03355675935745239]


fitness_maxl = [9.117834284204055e-13, 3.9615980465841678e-11, 9.4372981147492521e-11]
time_maxl = [0.15214191675186156, 0.04026339054107666, 0.03769617080688477]

fintess_nsteps = [2.6110851666361592, 885.88925568969694, 39.279918744711978]
time_nsteps = [0.037215137481689455, 0.004250586032867432, 0.017514777183532716]

fitness_modi = [1.3349758920682051, 18.196783986132637, 0.8653948686305728]
time_modi=[0.07368689775466919, 0.00827385187149048, 0.06883957386016845]


labels_clust = ["Same prameters","Modified parameters", "MaxLocIter 100", "NSteps=PopSize"]
labels_in = ["RMGS_all", "RMGS", "GGS"]

plt.rcParams['font.size'] = 12
fig, ax = plt.subplots()
ax.semilogy(time_same, fitness_same, 'o',  markersize=10)
ax.semilogy(time_modi, fitness_modi, 'ro', markersize=10)
ax.semilogy(time_maxl, fitness_maxl, 'yo', markersize=10)
ax.semilogy(time_nsteps, fintess_nsteps, 'go', markersize=10)
plt.legend(labels_clust)
for i, txt in enumerate(labels_in):
    ax.annotate(txt, (time_same[i],fitness_same[i]))
    ax.annotate(txt, (time_modi[i],fitness_modi[i]))
    ax.annotate(txt, (time_maxl[i],fitness_maxl[i]))
    ax.annotate(txt, (time_nsteps[i],fintess_nsteps[i]))
plt.draw()
plt.show()

# [3.6843911045038577e-11, 2.6159311516043271e-11, 0.19932895567739856]
# [0.03939963579177856, 0.03902599811553955, 0.0373329758644104]
# [3.4276472230117565e-11, 3.1061922786515307e-11, 3.25610341142305e-11]
# [0.03932163715362549, 0.039722990989685056, 0.037313199043273924]

# first results
# fitness_same = [2.13E-08, 311.387761409, 0.3747056007 ] 
# fitness_modi = [0, 1.62E-10, 0.0429593738]
# fitness_maxl = [0, 2.63E-11, 6.18E-12]
# fintess_nsteps = [0.064549041, 3756.55380419, 1.68E-10]
# time_same = [4.1913487911, 0.0644550323, 0.412115097]
# time_modi=[0, 0.282102108, 0.3806059361]
# time_maxl=[ 0, 0.5261609554, 0.5216841698]
# time_nsteps = [0.4833509922, 0.0189540386, 0.0940279961]


