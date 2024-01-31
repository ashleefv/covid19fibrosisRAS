import numpy as np
import math
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import ode
import random, pickle


############################################################
# User input for for RAS model
number_of_cases = 10
Nc = 2  # first loop to reach steady state without ACE2 (N1) and second loop with ACE2 (N2)
N1 = 1440  # min, simulation time before disease onset to load premorbid state
N2 = 14400  # min, simulation time after disease onset

#kT = 7e-4       # ANGII induced TGF-beta production rate
kT_all = [1e-3, 9e-4, 8e-4, 7e-4, 6e-4, 5e-4, 4e-4, 3e-4, 2e-4, 1e-4]       # ANGII induced TGF-beta production rate

pateint_type = 'hypertensive'
#pateint_type = 'normotensive'

#feedback_type = 'with'
feedback_type = 'without'
###############################################################

for Cases in range(number_of_cases):

    output = []
    result = []
    baseline = 0

    case = Cases + 1

    kT = kT_all[Cases]

    # hypertensive
    if (pateint_type == 'hypertensive'):
        x1 = 6e5        # Angiotensinogen  fmol/mL
        x2 = 110.0      # ANGI  fmol/mL
        x3 = 156.0      # ANGII  fmol/mL
        x4 = 92.0       # ANG(1-7)  fmol/mL
        x5 = 2.1e4      # ACE2-ANGII  fmol/mL
        x6 = 1.0        # ANGIV  fmol/mL
        x7 = 85.0       # ANGII-AT1R  fmol/mL
        x8 = 27.0       # ANGII-AT2R  fmol/mL

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
    x10 = 0.0           # Unconverted ANGII
    x11 = 0.0           # TGF-β  ng/mL
    x12 = 0.0           # Macrophage    (population number)
    x13 = 0.0           # Fibroblast    (population number)
    x14 = 0.0           # Collagen  μg



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

    fACE2 = np.array(pickle.load(open('ueACE2_'+str(1000)+'.p', 'rb')))
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
                if x[7]>=3:
                    beta = beta0 + (math.pow((x7/x[7]),delta)-1)
                else:
                    beta = beta0 + (math.pow((x7/3), delta) - 1)
            if (feedback_type == 'without'):
                beta = beta0


            fACE2 = np.array(pickle.load(open('ueACE2_'+str(1000)+'.p', 'rb')))
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
            y[10] = (kace2ang2on *fACE2[0]* x3 - kace2ang2on * xI * x[3])/(kace2ang2on *fACE2[0]* x3) -x[10]  # divide by kace2ang2on *fACE2[0]* x3 only to get the fraction for plotting

            # TGF-β
            y[11] = x[10] * kT *kace2ang2on *fACE2[0]* x3 + kMT * x[12] + kFT * x[13] - dT * x[11]

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


    result = np.concatenate((np.array(output[0]),np.array(output[1])))
    T = np.arange(-N1, N2, 0.1)


    #####################################################################################################
    # Saving files for plotting
    pickle.dump(result, open('resultT_KT'+str(Cases)+'.p', 'wb'))

