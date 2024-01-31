import numpy as np
import math
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
import random, pickle
from scipy.optimize import fsolve, least_squares
from scipy.optimize import curve_fit
from scipy.stats import qmc, norm, truncnorm, loguniform
from matplotlib.patches import Patch, Rectangle

ANGI_initial_samples = loguniform.rvs(2, 600, size=5000)
ANGII_initial_samples = loguniform.rvs(2, 1700, size=5000)
ANG17_initial_samples = loguniform.rvs(2, 150, size=5000)

ANGI_initial_value = []
ANGII_initial_value = []
ANG17_initial_value = []
sample_number = 1000
breakloop = 1
counter = 0
while breakloop>0:
    if ANGII_initial_samples[counter] >= 1.4*ANGI_initial_samples[counter] and ANGII_initial_samples[counter] >= 1.7*ANG17_initial_samples[counter] and ANGI_initial_samples[counter] >= 1.2*ANG17_initial_samples[counter]:
        ANGI_initial_value.append(ANGI_initial_samples[counter])
        ANGII_initial_value.append(ANGII_initial_samples[counter])
        ANG17_initial_value.append(ANG17_initial_samples[counter])
    counter = counter + 1
    if len(ANGII_initial_value) >= sample_number:
        print(counter)
        breakloop = -1

#print(ANGI_initial_value)
#print(len(ANGI_initial_value))
#input('')

'''
ANGI_initial_value = []
ANGII_initial_value = []
ANG17_initial_value = []

sampling_region = 10
sample_number = 10
sample_increment = 10
multiplier = 0.2
multiplier_increment = 0.2

sample_variation = 2
ANG17_higher_range = 2

for j in range(sample_variation):
    for i in range(sampling_region):
        breakloop = 1
        while breakloop>0:
            if j==0:
                mean = [69*multiplier, 166*multiplier, 9*multiplier]
            if j==1:
                mean = [69 * multiplier, 166 * multiplier, 9 * (multiplier+ANG17_higher_range)]
            std = [10,10,10]
            lhd = qmc.LatinHypercube(d=3).random(n=1)
            sample1 = norm(loc=mean, scale=std).ppf(lhd)
            if sample1[0,1] >= 1.4*sample1[0,0] and sample1[0,1] >= 1.7*sample1[0,2] and sample1[0,0] >= 1.2*sample1[0,2]:
                if sample1[0,0]>= 2 and sample1[0,0]<= 600 and sample1[0,1]>=2 and sample1[0,1]<=1700 and sample1[0,2]>=2 and sample1[0,2]<=200:
                    ANGI_initial_value.append(sample1[0, 0])
                    ANGII_initial_value.append(sample1[0,1])
                    ANG17_initial_value.append(sample1[0, 2])
                    #print(ANGII_initial_value)
                    #print(len(ANG17_initial_value))
            if len(ANGII_initial_value)>=sample_number:
                breakloop=-1
        multiplier = multiplier + multiplier_increment
        sample_number = sample_number + sample_increment
'''
ANGI = []
ANGII = []
ANG17 = []
ANGIIin = []
RAS_initial_value_list = []

# User input for for RAS model
number_of_cases = 1
Nc = 2  # first loop to reach steady state without ACE2 (N1) and second loop with ACE2 (N2)
N1 = 1440  # min, simulation time before disease onset to load premorbid state
N2 = 14400  # min, simulation time after disease onset
time_interval_from_ode_solution = 10        # depends on interval in ODE solver

kT = 5e-4  # ANGII induced TGF-beta production rate

for initial_value in range(len(ANGII_initial_value)):
    ANGI_value = ANGI_initial_value[initial_value]
    ANG_II_value = ANGII_initial_value[initial_value]
    ANG_17_value = ANG17_initial_value[initial_value]

    ############################################################
    pateint_type = 'hypertensive'
    #pateint_type = 'normotensive'

    #feedback_type = 'with'
    feedback_type = 'without'

    # saving the plot
    # Ones place (user-defined with cases): 1: HI, 2:NI, 3:HN, 4:NN
    # Tens place (fixed): 1: whole dynamics, 2: change from homeostasis with respect to ACE2
    counter1 = 11
    counter2 = 21

    ###############################################################

    for Cases in range(number_of_cases):

        output = []
        result = []
        baseline = 0

        case = Cases + 5

        # hypertensive
        if (pateint_type == 'hypertensive'):
            x1 = 6e5        # Angiotensinogen  fmol/mL
            x2 = ANGI_value      # ANGI  fmol/mL
            x3 = ANG_II_value      # ANGII  fmol/mL
            x4 = ANG_17_value       # ANG(1-7)  fmol/mL
            x5 = 2.1e4      # ACE2-ANGII  fmol/mL1
            x6 = 1.0*(1+(ANG_II_value-156)/156)        # ANGIV  fmol/mL
            x7 = 85.0*(1+(ANG_II_value-156)/156)        # ANGII-AT1R  fmol/mL
            #x7 = 50.4  # ANGII-AT1R  fmol/mL
            x8 = 27.0*(1+(ANG_II_value-156)/156)        # ANGII-AT2R  fmol/mL
            #x8 = 16.0  # ANGII-AT2R  fmol/mL

        # normotensive
        if (pateint_type == 'normotensive'):
            x1 = 6e5        # Angiotensinogen  fmol/mL
            x2 = 70         # ANGI  fmol/mL
            x3 = 28         # ANGII  fmol/mL
            x4 = 36         # ANG(1-7)  fmol/mL
            x5 = 2.1e4      # ACE2-ANGII  fmol/mL
            x6 = 1          # ANGIV  fmol/mL
            x7 = 15         # ANGII-AT1R  fmol/mL
            x8 = 5          # ANGII-AT2R  fmol/mL


        # Renin [x0] and MAS-ANG(1-7) [x9] are estimated using linear solving, see calculation after linear solving for details
        x10 = 0.0  # Unconverted ANGII
        x11 = 0.0  # TGF-β  ng/mL
        x12 = 0.0  # Macrophage    (population number)
        x13 = 0.0  # Fibroblast    (population number)
        x14 = 0.0  # Collagen  μg



        # parameters for RAS model
        hA = 600            # min
        hA1 = 0.5           # min
        hA2 = 0.5           # min
        hA17 = 0.5          # min
        hA4 = 0.5           # min
        hAT1 = 12           # min
        hAT2 = 12           # min
        hmas = 12           # min
        hR = 12             # min
        cR = 20             # 1/min
        delta = 0.8         # dimensionless parameter to account for the effect of downstream feedback in RAS network


        # parameters for Immune model
        dM = 0.6 / (24 * 60)
        kF = 0.924 / (24 * 60)
        dF = 0.12 / (24 * 60)
        dTB = 15 / (24 * 60)
        kFC = 2.52e-7           # μg/min

        # TGF-beta
        dT = 1.04e-2
        # parameters from Jin and ceresa. Their figures are plotted for microliter. To convert pg to ng and uL to mL, units will be similar.

        kMT = 0.07 / (24 * 60)      # TGF-β secretion from macrophages (ng/cell/day in mL (converted from pg/cell/day))
        kFT = 0.004 / (24 * 60)     # TGF-β secretion from fibroblasts

        fACE2 = np.array(pickle.load(open('ueACE2_'+str(case*200)+'.p', 'rb')))
        print(fACE2[0])


        ##################################################################################################################
        # Linear Solver to solve following set of equations using steady state initial values of RAS peptides

        # beta0 - math.log(2) / hR * x0 = 0
        # kA - cR * x0 - math.log(2) / hA * x1 = 0
        # cR * x0 - cA * ACE_0 * x2 - cN * x2 - math.log(2) / hA1 * x2 = 0
        # cA * ACE_0 * x2 - cAT1 * x3 - cAT2 * x3 - cAPA * x3 - kace2ang2on * ACE2_0 * x3 - math.log(2) / hA2 * x3 = 0
        # cN * x2 + kace2 * x5 - cmas * x4 - math.log(2) / hA17 * x4 = 0
        # kace2ang2on * ACE2_0 * x3 - kace2 * x5 = 0
        # cAPA * x3 - math.log(2) / hA4 * x6 = 0
        # cAT1 * x3 - math.log(2) / hAT1 * x7 = 0
        # cAT2 * x3 - math.log(2) / hAT2 * x8 = 0
        # cmas * x4 - math.log(2) / hmas * x9 = 0

        # 10 unknown parameters are: x0, beta0, kA, cA, kace2ang2on, kace2, cAPA, cAT1, cAT2, x9
        # model assumptions for additional parameters: cmas = cAT2,  cN = 0


        a = np.array([[-math.log(2) / hR, 1, 0, 0, 0, 0, 0, 0, 0, 0], [-cR, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                      [cR, 0, 0, -x2* fACE2[0], 0, 0, 0, 0, 0, 0], [0, 0, 0, x2* fACE2[0],  -x3 * fACE2[0], 0, -x3, -x3, -x3, 0],
                      [0, 0, 0, 0, 0, x5, 0, 0, -x4, 0], [0, 0, 0, 0, x3 * fACE2[0], -x5, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, x3, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, x3, 0, 0],
                      [0, 0, 0, 0, 0, 0, 0, 0, x3, 0], [0, 0, 0, 0, 0, 0, 0, 0, x4, -math.log(2) / hmas]])
        b = np.array([0,math.log(2)/hA * x1, math.log(2)/hA1 * x2, math.log(2)/hA2* x3 , math.log(2)/hA17 * x4, 0, math.log(2)/hA4 * x6, math.log(2)/hAT1 * x7, math.log(2)/hAT2 * x8,0])

        solution = np.linalg.solve(a, b)
        print(solution)

        x0 = solution[0]                # Renin  fmol/mL
        beta0 = solution[1]             # fmol/(mL min)
        kA = solution[2]                # fmol/(mL min)
        cA = solution[3]                # 1/(min (# of ACE receptors)), calibration with virtual lung tissue
        kace2ang2on = solution[4]       # 1/(min (# of ACE2 receptors)), calibration with virtual lung tissue
        kace2 = solution[5]             # 1/min
        cAPA = solution[6]              # 1/min
        cAT1 = solution[7]              # 1/min
        cAT2 = solution[8]              # 1/min
        x9 = solution[9]                # MAS-ANG(1-7)  fmol/mL

        # model assumptions
        cmas = cAT2
        cN = 0

        ###########################################################################################


        for i in range(Nc):
            def diff(x, T, baseline,kT):

                if (feedback_type == 'with'):
                    if x[7]>3:
                        beta = beta0 + (math.pow((x7/x[7]),delta)-1)
                    else:
                        beta = beta0 + (math.pow((x7/3), delta) - 1)
                if (feedback_type == 'without'):
                    beta = beta0


                fACE2 = np.array(pickle.load(open('ueACE2_'+str(case*200)+'.p', 'rb')))
                tACE2 = np.array(pickle.load(open('ueACE2t.p', 'rb')))


                if baseline ==0:
                    xI = fACE2[0]
                if baseline == 1:
                    xI = np.interp((T/(24*60)), tACE2, fACE2)



                ##########################################################################################################
                # interpolation and extrapolation from experimental data

                # Collected experimental data from literature for macrophage recruitment
                M_data = [30.48, 61.34, 125.78, 88.99, 81.28, 56.42, 46.70]
                T_beta_data1 = [0.05/1000, 0.1/1000, 0.5/1000, 1/1000, 10/1000, 100/1000, 1000/1000]

                # Collected experimental data from literature for fibroblasts recruitment
                Fg_data = [5.35, 8.41, 8.58, 13.58, 28.83]
                T_beta_data2 = [0, 0.01, 0.1, 1, 10]

                # Considering the difference from interpolated value to baseline value
                MTB = np.interp(x[11], T_beta_data1, M_data) - M_data[0]
                FgTB = np.interp(x[11], T_beta_data2, Fg_data) - Fg_data[0]

                ##########################################################################################################


                # Differential equations
                y = [0.0 for i in range(len(x))]

                # Renin
                y[0] = beta - math.log(2) / hR * x[0]

                # Angiotensinogen
                y[1] = kA - cR * (xI/fACE2[0]) * x[0] - math.log(2) / hA * x[1]

                # ANGI
                y[2] = cR * (xI/fACE2[0]) * x[0] - cA * xI * x[2] - cN * x[2] - math.log(2) / hA1 * x[2]

                # ANGII
                y[3] = cA * xI * x[2] - cAT1 * x[3] - cAT2 * x[3] - cAPA * x[3] - kace2ang2on * xI * x[3] - math.log(2) / hA2 * x[3]

                # ANG(1-7)
                y[4] = cN * x[2] + kace2 * x[5] - cmas * x[4] - math.log(2) / hA17 * x[4]

                # ACE2-ANGII
                y[5] = kace2ang2on * xI * x[3] - kace2 * x[5]

                # ANGIV
                y[6] = cAPA * x[3] - math.log(2) / hA4 * x[6]

                # ANGII-AT1R
                y[7] = cAT1 * x[3] - math.log(2) / hAT1 * x[7]

                # ANGII-AT2R
                y[8] = cAT2 * x[3] - math.log(2) / hAT2 * x[8]

                # MAS-ANG(1-7)
                y[9] = cmas * x[4] - math.log(2) / hmas * x[9]

                # Unconverted ANGII
                y[10] = (kace2ang2on * fACE2[0] * x3 - kace2ang2on * xI * x[3]) / (kace2ang2on * fACE2[0] * x3) - x[10]  # divide by kace2ang2on *fACE2[0]* x3 only to get the fraction for plotting

                # TGF-β
                y[11] = x[10] * kT * kace2ang2on * fACE2[0] * x3 + kMT * x[12] + kFT * x[13] - dT * x[11]

                # Macrophage
                y[12] = MTB / 90 - dM * x[12]

                # Fibroblast
                y[13] = FgTB / (48 * 60) - dF * x[13]

                # Collagen
                y[14] = kFC * ((0.942 * x[11]) / (0.174 + x[11])) * x[13]

                return y


            if i == 0:
                x = (x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14)
                T = np.arange(0.0, N1, 0.1)
                result = np.array(odeint(diff, x, T, args=(0, kT)))
                output.append(result)

            else:
                x = (result[-1, 0], result[-1, 1], result[-1, 2], result[-1, 3], result[-1, 4], result[-1, 5], result[-1, 6], result[-1, 7], result[-1, 8], result[-1, 9], result[-1, 10], result[-1, 11], result[-1, 12], result[-1, 13], result[-1, 14])
                T = np.arange(0.0, N2, 0.1)
                result = np.array(odeint(diff, x, T, args=(1, kT)))
                output.append(result)

        result = np.concatenate((np.array(output[0]), np.array(output[1])))
        T = np.arange(-N1, N2, 0.1)

        #####################################################################################################
        # Saving files for plotting
        pickle.dump(result, open('resultT' + str(Cases) + '.p', 'wb'))

    #######################################################################################################
    # plotting section

    N3 = N1*time_interval_from_ode_solution
    tc = np.arange(0, N2, 0.1)      # from disease onset
    tc1 = np.arange(-N1, N2, 0.1)   # from preload condition


    resultC1 = [[] for _ in range(number_of_cases)]      # dysregulated ACE2
    ueACE2 = [[] for _ in range(number_of_cases)]        # cell surface ACE2 dynamics

    for i in range(number_of_cases):
        result = np.array(pickle.load(open('resultT'+str(i)+'.p', 'rb')))
        resultC1[i].append(result)

    resultC1 = np.array(resultC1)

    for i in range(number_of_cases):
        result = np.array(pickle.load(open('ueACE2_'+str((i+1)*200)+'.p', 'rb')))
        ueACE2[i].append(result)

    ueACE2 = np.array(ueACE2)
    t1 = np.array(pickle.load(open('ueACE2t.p', 'rb')))

    ANGI.extend([resultC1[0, 0, :, 2]])
    ANGII.extend([resultC1[0, 0, :, 3]])
    ANG17.extend([resultC1[0, 0, :, 4]])
    ANGIIin.extend([resultC1[0, 0, :, 10]*100])

# Saving files for plotting
pickle.dump(ANGI, open('ANGI.p', 'wb'))
pickle.dump(ANGII, open('ANGII.p', 'wb'))
pickle.dump(ANG17, open('ANG17.p', 'wb'))
pickle.dump(ANGIIin, open('ANGIIin.p', 'wb'))

plt.rcParams.update({'font.size': 25})
N3 = N1*time_interval_from_ode_solution
tc = np.arange(0, N2, 0.1)      # from disease onset
tc1 = np.arange(-N1, N2, 0.1)   # from preload condition


# ANGII
fig, ax = plt.subplots()
for i in range(len(ANGII)):
    Value_ANGI_array = np.array(ANGI[i])
    last_value_ANGI = Value_ANGI_array[-1]
    Value_ANGII_array = np.array(ANGII[i])
    last_value_ANGII = Value_ANGII_array[-1]
    Value_ANG17_array = np.array(ANG17[i])
    last_value_ANG17 = Value_ANG17_array[-1]
    if last_value_ANGII >= 2.1 and last_value_ANGII <= 357 and last_value_ANG17 >= 3 and last_value_ANG17 <= 7.3 and last_value_ANGI >= 2.7 and last_value_ANGI <= 188:
        ax.plot(tc1 / (24 * 60), np.array(ANGII[i]), linewidth=2)
        RAS_initial_value_list.extend([[Value_ANGI_array[0], Value_ANGII_array[0], Value_ANG17_array[0]]])

print(RAS_initial_value_list)
print(len(RAS_initial_value_list))
# experimental data from Kutz
positions=[-1, 10]
boxes = [
    {
        'whislo': 10,    # Bottom whisker position
        'q1'    : 49,    # First quartile (25th percentile)
        'med'   : 166,    # Median         (50th percentile)
        'q3'    : 367,    # Third quartile (75th percentile)
        'whishi': 1687,    # Top whisker position
        'fliers': [2.4]        # Outliers
    },
    {
        'whislo': 2.1,    # Bottom whisker position
        'q1'    : 17.3,    # First quartile (25th percentile)
        'med'   : 38,    # Median         (50th percentile)
        'q3'    : 114,    # Third quartile (75th percentile)
        'whishi': 357,    # Top whisker position
        'fliers': []        # Outliers
    },

]
ax.bxp(boxes, positions, manage_ticks = False, widths=1)

plt.legend()
old_handles, labels = ax.get_legend_handles_labels()
legend_elements = [Patch(facecolor='white', edgecolor='k', label='Kutz et al.')]

plt.legend(handles=old_handles + legend_elements, loc='upper right', fontsize="15")

#plt.legend(loc='lower left', prop={"size": 12})
#plt.title('Kutz et al. ')
plt.xlabel('Time (day)')
plt.ylabel('ANGII')
plt.yscale("log")
#plt.ylim([-105, 5])
plt.savefig("ANGII_values_kutz.png", dpi=300, bbox_inches='tight')
plt.show()

# ANGII full list
fig0, ax0 = plt.subplots()
for i in range(len(ANGII)):
    ax0.plot(tc1 / (24 * 60), np.array(ANGII[i]), linewidth=2)

ax0.bxp(boxes, positions, manage_ticks = False, widths=1)

plt.legend()
old_handles, labels = ax0.get_legend_handles_labels()
legend_elements = [Patch(facecolor='white', edgecolor='k', label='Kutz et al.')]

plt.legend(handles=old_handles + legend_elements, loc='upper right', fontsize="15")

#plt.legend(loc='lower left', prop={"size": 12})
#plt.title('Kutz et al. ')
plt.xlabel('Time (day)')
plt.ylabel('ANGII')
plt.yscale("log")
#plt.ylim([-105, 5])
plt.savefig("ANGII_values_kutz_full.png", dpi=300, bbox_inches='tight')
plt.show()

# ANG(1-7)
fig1, ax1 = plt.subplots()
for i in range(len(ANG17)):
    Value_ANGI_array = np.array(ANGI[i])
    last_value_ANGI = Value_ANGI_array[-1]
    Value_ANGII_array = np.array(ANGII[i])
    last_value_ANGII = Value_ANGII_array[-1]
    Value_ANG17_array = np.array(ANG17[i])
    last_value_ANG17 = Value_ANG17_array[-1]
    if last_value_ANGII >= 2.1 and last_value_ANGII <= 357 and last_value_ANG17 >= 3 and last_value_ANG17 <= 7.3 and last_value_ANGI >= 2.7 and last_value_ANGI <= 188:
        ax1.plot(tc1 / (24 * 60), np.array(ANG17[i]), linewidth=2)

# experimental data from Kutz
positions=[-1, 10]
boxes = [
    {
        'whislo': 3,    # Bottom whisker position
        'q1'    : 3,    # First quartile (25th percentile)
        'med'   : 7.6,    # Median         (50th percentile)
        'q3'    : 9,    # Third quartile (75th percentile)
        'whishi': 38.1,    # Top whisker position
        'fliers': [58, 111, 166]        # Outliers
    },
    {
        'whislo': 3,    # Bottom whisker position
        'q1'    : 3,    # First quartile (25th percentile)
        'med'   : 3,    # Median         (50th percentile)
        'q3'    : 3,    # Third quartile (75th percentile)
        'whishi': 3,    # Top whisker position
        'fliers': [3.3, 4.9, 5.4, 7, 7.3]        # Outliers
    },

]
ax1.bxp(boxes, positions, manage_ticks = False, widths=1)

plt.legend()
old_handles, labels = ax1.get_legend_handles_labels()
legend_elements = [Patch(facecolor='white', edgecolor='k', label='Kutz et al.')]

plt.legend(handles=old_handles + legend_elements, loc='upper right', fontsize="15")

#plt.legend(loc='lower left', prop={"size":12})
#plt.title('Kutz et al. ')
plt.xlabel('Time (day)')
plt.ylabel('ANG1-7')
plt.yscale("log")
#plt.ylim([-105,5])
plt.savefig("ANG1-7_values_kutz.png", dpi = 300, bbox_inches='tight')
plt.show()


# ANG(1-7) full list
fig11, ax11 = plt.subplots()
for i in range(len(ANG17)):
    ax11.plot(tc1 / (24 * 60), np.array(ANG17[i]), linewidth=2)

ax11.bxp(boxes, positions, manage_ticks = False, widths=1)

plt.legend()
old_handles, labels = ax11.get_legend_handles_labels()
legend_elements = [Patch(facecolor='white', edgecolor='k', label='Kutz et al.')]

plt.legend(handles=old_handles + legend_elements, loc='upper right', fontsize="15")

#plt.legend(loc='lower left', prop={"size":12})
#plt.title('Kutz et al. ')
plt.xlabel('Time (day)')
plt.ylabel('ANG1-7')
plt.yscale("log")
#plt.ylim([-105,5])
plt.savefig("ANG1-7_values_kutz_full.png", dpi = 300, bbox_inches='tight')
plt.show()


# ANGI
fig2, ax2 = plt.subplots()
for i in range(len(ANGI)):
    Value_ANGI_array = np.array(ANGI[i])
    last_value_ANGI = Value_ANGI_array[-1]
    Value_ANGII_array = np.array(ANGII[i])
    last_value_ANGII = Value_ANGII_array[-1]
    Value_ANG17_array = np.array(ANG17[i])
    last_value_ANG17 = Value_ANG17_array[-1]
    if last_value_ANGII >= 2.1 and last_value_ANGII <= 357 and last_value_ANG17 >= 3 and last_value_ANG17 <= 7.3 and last_value_ANGI >= 2.7 and last_value_ANGI <= 188:
        ax2.plot(tc1 / (24 * 60), np.array(ANGI[i]), linewidth=2)

# experimental data from Kutz
positions=[-1, 10]
boxes = [
    {
        'whislo': 3,    # Bottom whisker position
        'q1'    : 8,    # First quartile (25th percentile)
        'med'   : 69,    # Median         (50th percentile)
        'q3'    : 149,    # Third quartile (75th percentile)
        'whishi': 595,    # Top whisker position
        'fliers': []        # Outliers
    },
    {
        'whislo': 2.7,    # Bottom whisker position
        'q1'    : 5,    # First quartile (25th percentile)
        'med'   : 15,    # Median         (50th percentile)
        'q3'    : 53,    # Third quartile (75th percentile)
        'whishi': 188,    # Top whisker position
        'fliers': []        # Outliers
    },

]
ax2.bxp(boxes, positions, manage_ticks = False, widths=1)

plt.legend()
old_handles, labels = ax2.get_legend_handles_labels()
legend_elements = [Patch(facecolor='white', edgecolor='k', label='Kutz et al.')]

plt.legend(handles=old_handles + legend_elements, loc='upper right', fontsize="15")

#plt.legend(loc='lower left', prop={"size": 12})
#plt.title('Kutz et al. ')
plt.xlabel('Time (day)')
plt.ylabel('ANGI')
plt.yscale("log")
#plt.ylim([-105, 5])
plt.savefig("ANGI_values_kutz.png", dpi=300, bbox_inches='tight')
plt.show()

# ANGI full list
fig21, ax21 = plt.subplots()
for i in range(len(ANGI)):
    ax21.plot(tc1 / (24 * 60), np.array(ANGI[i]), linewidth=2)

ax21.bxp(boxes, positions, manage_ticks = False, widths=1)

plt.legend()
old_handles, labels = ax21.get_legend_handles_labels()
legend_elements = [Patch(facecolor='white', edgecolor='k', label='Kutz et al.')]

plt.legend(handles=old_handles + legend_elements, loc='upper right', fontsize="15")

#plt.legend(loc='lower left', prop={"size": 12})
#plt.title('Kutz et al. ')
plt.xlabel('Time (day)')
plt.ylabel('ANGI')
plt.yscale("log")
#plt.ylim([-105, 5])
plt.savefig("ANGI_values_kutz_full.png", dpi=300, bbox_inches='tight')
plt.show()

# Unconverted ANGII

'''
# Baseline experimental range
Baseline_day = [-1]
Baseline_mean = [77]
Baseline_error = [52]
plt.errorbar(Baseline_day, Baseline_mean, yerr=Baseline_error, color = 'green', fmt='o',capsize=4, label='Baseline')

# Non Critically ill
Non_critical_day = [7]
Non_critical_mean = [127.4]
Non_critical_error = [18.2]
plt.errorbar(Non_critical_day, Non_critical_mean, yerr=Non_critical_error, color = 'red', fmt='.',capsize=4, label='Non-critical')

# Critically ill
Critical_day = [7]
Critical_mean = [147.4]
Critical_error = [22.1]
plt.errorbar(Critical_day, Critical_mean, yerr=Critical_error, color = 'red', fmt='o',capsize=4, label='Critical')
'''

# Unconverted ANGII
fig3, ax3 = plt.subplots()
normal_mean = (25+129)/2
upper_th = ((169.5-normal_mean)/normal_mean)*100
lower_th = ((131.5-normal_mean)/normal_mean)*100
mean_range = ((147.4-normal_mean)/normal_mean)*100


#upper_th = ((169.5-127.4)/127.4)*100
#lower_th = ((131.5-127.4)/127.4)*100
#mean_range = ((147.4-127.4)/127.4)*100
# Fill in area under the curve and the horizontal lines
ax3.fill_between(x=tc1 / (24 * 60), y1=upper_th, y2=lower_th, color='blue',  interpolate=True, alpha=.05)

# Horizontal lines
#plt.axhline(upper_th, linestyle='-.', color='black', label="maximum")
ax3.axhline(mean_range, linestyle='-.', color='red', label='Liu et al.', linewidth=2)
#plt.axhline(lower_th, linestyle='-.', color='orange', label='minimum')

for i in range(len(ANGIIin)):
    ax3.plot(tc1 / (24 * 60), np.array(ANGIIin[i]), linewidth=2)
plt.legend(loc='lower right', prop={"size":15})
#old_handles, labels = ax.get_legend_handles_labels()
#legend_elements = [Patch(facecolor='lightskyblue', edgecolor='k', label='Liu et al.')]
#plt.legend(handles=old_handles + legend_elements, loc='lower right', fontsize="15")
#plt.title('Liu et al. ')
plt.xlabel('Time (day)')
plt.ylabel('% Unconverted ANGII')
#plt.ylim([-2.5,55])
#plt.yscale("log")
plt.savefig("ANGII_change_Liu.png", dpi = 300, bbox_inches='tight')
plt.show()



