# -*- coding: utf-8 -*-
"""
Created on Thu May  6 17:53:32 2021

@author: d338c921
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats

abundances = pd.read_csv('Brewer_Fischer_abundance_ratios.csv')
C_O = abundances.values[:, 7].T
Mg_Si = abundances.values[:, 8].T
Np = abundances.values[:, 9].T

C_O = C_O.astype('float')
Mg_Si = Mg_Si.astype('float')
Np = Np.astype('float')

#gaussian fits
mu_C_O, sig_C_O = stats.norm.fit(C_O)
mu_Mg_Si, sig_Mg_Si = stats.norm.fit(Mg_Si)
n_C_O, x_C_O  = np.histogram(C_O, bins = 100)
n_Mg_Si, x_Mg_Si = np.histogram(Mg_Si, bins = 100)
y_C_O = stats.norm.pdf(x_C_O, mu_C_O, sig_C_O)
y_Mg_Si = stats.norm.pdf(x_Mg_Si, mu_Mg_Si, sig_Mg_Si)

#repeat the process above for only exoplanet hosts
host_C_O = []
host_Mg_Si = []
for i in range(len(Np)):
    if (np.isnan(Np[i]) == False):
        host_C_O = np.append(host_C_O, C_O[i])
        host_Mg_Si = np.append(host_Mg_Si, Mg_Si[i])
        
mu_C_O_host, sig_C_O_host = stats.norm.fit(host_C_O)
mu_Mg_Si_host, sig_Mg_Si_host = stats.norm.fit(host_Mg_Si)
n, x_C_O_host = np.histogram(host_C_O, bins = 100)
n, x_Mg_Si_host = np.histogram(host_Mg_Si, bins = 100)
y_C_O_host = stats.norm.pdf(x_C_O_host, mu_C_O_host, sig_C_O_host)
y_Mg_Si_host = stats.norm.pdf(x_Mg_Si_host, mu_Mg_Si_host, sig_Mg_Si_host)
        
        
plt.figure()
plt.title("Distribution of C/O Abundance Ratios")
plt.xlabel("C/O Ratio")
plt.ylabel("Number of Stars")
plt.hist(C_O, bins = 100, color = "dodgerblue", label = "All Stars")
plt.hist(host_C_O, bins = 100, color = "r", label = "Host Stars")
#plot solar values for C/O and Mg/Si
plt.axvline(x = 0.54, color = "darkorange", linestyle = "--", label = "Solar Value: 0.54") #From Brewer & Fischer 2016
plt.axvline(x = mu_C_O, color = "dodgerblue", linestyle = "--", label = "Mean All Star Ratio: " + str(np.round(mu_C_O, 2)))
plt.axvline(x = mu_C_O_host, color = "r", linestyle = "--", label = "Mean Host Ratio: " + str(np.round(mu_C_O_host, 2)))
plt.legend()

plt.figure()
plt.title("Distribution of Mg/Si Abundance Ratios")
plt.xlabel("Mg/Si Ratio")
plt.ylabel("Number of Stars")
plt.hist(Mg_Si, bins = 100, color = "dodgerblue", label = "All Stars")
plt.hist(host_Mg_Si, bins = 100, color = "r", label = "Host Stars")
#plot solar values for Mg/Si
plt.axvline(x = 1.05, color = "darkorange", linestyle = "--", label = "Solar Value: 1.05") #From Brewer & Fischer 2016
plt.axvline(x = mu_Mg_Si, color = "dodgerblue", linestyle = "--", label = "Mean All Star Ratio: " + str(np.round(mu_Mg_Si, 2)))
plt.axvline(x = mu_Mg_Si_host, color = "r", linestyle = "--", label = "Mean Host Ratio: " + str(np.round(mu_Mg_Si_host, 2)))
plt.legend()

plt.figure()
plt.title("Distribution of C/O Abundance Ratios: Models")
plt.xlabel("C/O Ratio")
plt.ylabel("Probability Distribution Function")
plt.plot(x_C_O, y_C_O, color = "dodgerblue", label = "All Star Model")
plt.plot(x_C_O_host, y_C_O_host, color = "r", label = "Host Star Model")
plt.axvline(x = 0.54, color = "darkorange", linestyle = "--", label = "Solar Value: 0.54") #From Brewer & Fischer 2016
plt.axvline(x = mu_C_O, color = "dodgerblue", linestyle = "--", label = "Mean All Star Ratio: " + str(np.round(mu_C_O, 2)))
plt.axvline(x = mu_C_O_host, color = "r", linestyle = "--", label = "Mean Host Ratio: " + str(np.round(mu_C_O_host, 2)))
plt.legend()

plt.figure()
plt.title("Distribution of Mg/Si Abundance Ratios: Models")
plt.xlabel("Mg/Si Ratio")
plt.ylabel("Probability Distribution Function")
plt.plot(x_Mg_Si, y_Mg_Si, color = "dodgerblue", label = "All Star Model")
plt.plot(x_Mg_Si_host, y_Mg_Si_host, color = "r", label = "Host Star Model")
plt.axvline(x = 1.05, color = "darkorange", linestyle = "--", label = "Solar Value: 1.05") #From Brewer & Fischer 2016
plt.axvline(x = mu_Mg_Si, color = "dodgerblue", linestyle = "--", label = "Mean All Star Ratio: " + str(np.round(mu_Mg_Si, 2)))
plt.axvline(x = mu_Mg_Si_host, color = "r", linestyle = "--", label = "Mean Host Ratio: " + str(np.round(mu_Mg_Si_host, 2)))
plt.legend()

###
#plot horizontal lines located at each of the means
###

#Use cumulative distribution functions to calculte the probability that a star has sub- or super-solar values
# plt.figure()
# plt.title("Distribution of C/O Abundance Ratios: Models")
# plt.xlabel("x")
# plt.ylabel("Cumulative Distribution Function: F(x) = P(X <= x)")
# plt.plot(x_C_O, stats.norm.cdf(x_C_O, mu_C_O, sig_C_O), label = "All Star Model")
# plt.plot(x_C_O_host, stats.norm.cdf(x_C_O_host, mu_C_O_host, sig_C_O_host), label = "Host Star Model")

# plt.figure()
# plt.title("Distribution of Mg/Si Abundance Ratios: Models")
# plt.xlabel("x")
# plt.ylabel("Cumulative Distribution Function: F(x) = P(X <= x)")
# plt.plot(x_Mg_Si, stats.norm.cdf(x_Mg_Si, mu_Mg_Si, sig_Mg_Si), label = "All Star Model")
# plt.plot(x_Mg_Si_host, stats.norm.cdf(x_Mg_Si_host, mu_Mg_Si_host, sig_Mg_Si_host), label = "Host Star Model")

#Don't really need plots of the CDFs, just need them for calculations
prob_subsol_all_CO = stats.norm.cdf(0.54, mu_C_O, sig_C_O)
prob_supersol_all_CO = 1 - prob_subsol_all_CO

prob_subsol_host_CO = stats.norm.cdf(0.54, mu_C_O_host, sig_C_O_host)
prob_supersol_host_CO = 1 - prob_subsol_host_CO

prob_subsol_all_MgSi = stats.norm.cdf(1.05, mu_Mg_Si, sig_Mg_Si)
prob_supersol_all_MgSi = 1 - prob_subsol_all_MgSi

prob_subsol_host_MgSi = stats.norm.cdf(1.05, mu_Mg_Si_host, sig_Mg_Si_host)
prob_supersol_host_MgSi = 1 - prob_subsol_host_MgSi

print("The probability of finding a star with a sub-solar C/O ratio (among all stars) is: " + str(np.round(prob_subsol_all_CO, 2) * 100) + "%\n")
print("The probability of finding a star with a super-solar C/O ratio (among all stars) is: " + str(np.round(prob_supersol_all_CO, 2) * 100) + "%\n")
print("The probability of finding a star with a sub-solar C/O ratio (among host stars) is: " + str(np.round(prob_subsol_host_CO, 2) * 100) + "%\n")
print("The probability of finding a star with a super-solar C/O ratio (among host stars) is: " + str(np.round(prob_supersol_host_CO, 2) * 100) + "%\n")

print("The probability of finding a star with a sub-solar Mg/Si ratio (among all stars) is: " + str(np.round(prob_subsol_all_MgSi, 2) * 100) + "%\n")
print("The probability of finding a star with a super-solar Mg/Si ratio (among all stars) is: " + str(np.round(prob_supersol_all_MgSi, 2) * 100) + "%\n")
print("The probability of finding a star with a sub-solar Mg/Si ratio (among host stars) is: " + str(np.round(prob_subsol_host_MgSi, 2) * 100) + "%\n")
print("The probability of finding a star with a super-solar Mg/Si ratio (among host stars) is: " + str(np.round(prob_supersol_host_MgSi, 2) * 100) + "%\n")

P_CO = stats.norm.cdf(0.8, mu_C_O, sig_C_O) - stats.norm.cdf(0.5, mu_C_O, sig_C_O)

P_MgSi = stats.norm.cdf(2.0, mu_Mg_Si, sig_Mg_Si) - stats.norm.cdf(1.0, mu_Mg_Si, sig_Mg_Si)

print("Probability of finding a star with Earth-like C/O: " + str(np.round(P_CO * 100, 3)))
print("Probability of finding a star with Earth-like Mg/Si: " + str(np.round(P_MgSi * 100, 3)))
print("Probability of finding a star with Earth-like C/O AND Mg/Si: " + str(np.round(P_CO * P_MgSi * 100, 3)))

P_CO_host = stats.norm.cdf(0.8, mu_C_O_host, sig_C_O_host) - stats.norm.cdf(0.5, mu_C_O_host, sig_C_O_host)

P_MgSi_host = stats.norm.cdf(2.0, mu_Mg_Si_host, sig_Mg_Si_host) - stats.norm.cdf(1.0, mu_Mg_Si_host, sig_Mg_Si_host)

print("Among host stars:")
print("Probability of finding a star with Earth-like C/O: " + str(np.round(P_CO_host * 100, 3)))
print("Probability of finding a star with Earth-like Mg/Si: " + str(np.round(P_MgSi_host * 100, 3)))
print("Probability of finding a star with Earth-like C/O AND Mg/Si: " + str(np.round(P_CO_host * P_MgSi_host * 100, 3)))


