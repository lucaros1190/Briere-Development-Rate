# Python script to check the values of T_opt in Briére function
# Suggestion provided by Reviewer 1
# Created by Luca Rossini on 26 April 2022
# E-mail: luca.rossini@unitus.it
# Last update 26 April 2022

from math import *

# Definition of the function to calculate T_opt in Briére function

def ToptBriere(a, T_L, T_M, m):
    T_opt = T_opt = (2 * m * T_M + T_L * (m + 1) + sqrt(4 * m * m * T_M * T_M + T_L * T_L * (m + 1) * (m + 1) - 4 * m * m * T_M * T_L))  / (4 * m + 2)
    return T_opt

# Definition of the function to calculate the error associated with T_opt in Briére function

def Err_ToptBriere(a, ea, T_L, eT_L, T_M, eT_M, m, em):

    A = pow(eT_L, 2) * pow( (m + 1)/(4*m + 2) + (1/(4*m + 2)) * ((T_L * pow((m + 1), 2) - 2*m * T_M)/(sqrt(4 * m * m * T_M * T_M + T_L * T_L * pow(m+1, 2) - 4*m * T_M * T_L))) ,2);
                
    B = pow(eT_L, 2) * pow(m/(2*m +1) + (1/(4*m + 2) * (4*m*m * T_M - 2*m * T_L)/(sqrt(4*m*m * T_M*T_M + T_L*T_L * pow((m+1), 2) - 4 * m * m * T_L * T_M))) ,2);
                
    C = pow(em, 2) * pow(2*T_M * ((1/(4*m + 2)) - ((4*m)/(pow((4*m + 2) ,2))) + (T_L/((pow((4*m + 2) ,2)))) * ((4*T_M*T_M * (m - 2*m*m) - T_L*T_L * (m+1) * (4*m + 3) + T_M * T_M * (8*m - 2) )/( (pow((4*m + 2) ,2)) * sqrt(4*m*m * T_M*T_M + T_L*T_L * pow((m + 1), 2) - 4*m*T_L * T_M) )) ) ,2);
    
    eT_opt = sqrt(A + B + C)
    
    return eT_opt


# Calculation of optimal temperatures and associated errors for all the isolates involved into the study

    # Isolate P-1

a_P1 = 0.0063342
T_L_P1 = -2.2291
T_M_P1 = 29.7919
m_P1 = 1.9342086
ea_P1 = 0.002292
eT_L_P1 = 4.4065
eT_M_P1 = 1.0718
em_P1 = 0.7250

T_optP1 = ToptBriere(a_P1, T_L_P1, T_M_P1, m_P1)
Err_ToptP1 = Err_ToptBriere(a_P1, ea_P1, T_L_P1, eT_L_P1, T_M_P1, eT_M_P1, m_P1, em_P1)

print('\n  Isolate P-1:\n')
print('  T_opt =', T_optP1, '+\-', Err_ToptP1, '\n')


    # Isolate T-16

a_T16 = 0.00343
T_L_T16 = -8.4479
T_M_T16 = 31.3087
m_T16 = 1.4502
ea_T16 = 0.001
eT_L_T16 = 5.6696
eT_M_T16 = 0.9765
em_T16 = 0.4333

T_optT16 = ToptBriere(a_T16, T_L_T16, T_M_T16, m_T16)
Err_ToptT16 = Err_ToptBriere(a_T16, ea_T16, T_L_T16, eT_L_T16, T_M_T16, eT_M_T16, m_T16, em_T16)

print('  Isolate T-16:\n')
print('  T_opt =', T_optT16, '+\-', Err_ToptT16, '\n')


    # Isolate F-1C

a_F1C = 0.002264
T_L_F1C = 0
T_M_F1C = 31.9196
m_F1C = 1.0496
ea_F1C = 0.001068
eT_L_F1C = 0
eT_M_F1C = 0.9894
em_F1C = 0.1663

T_optF1C = ToptBriere(a_F1C, T_L_F1C, T_M_F1C, m_F1C)
Err_ToptF1C = Err_ToptBriere(a_F1C, ea_F1C, T_L_F1C, eT_L_F1C, T_M_F1C, eT_M_F1C, m_F1C, em_F1C)

print('  Isolate F-1C:\n')
print('  T_opt =', T_optF1C, '+\-', Err_ToptF1C, '\n')


    # Isolate TG-1

a_TG1 = 0.0032938
T_L_TG1 = -0.5391
T_M_TG1 = 31.117
m_TG1 = 1.1706
ea_TG1 = 0.0009145
eT_L_TG1 = 1.2897
eT_M_TG1 = 0.492
em_TG1 = 0.1625

T_optTG1 = ToptBriere(a_TG1, T_L_TG1, T_M_TG1, m_TG1)
Err_ToptTG1 = Err_ToptBriere(a_TG1, ea_TG1, T_L_TG1, eT_L_TG1, T_M_TG1, eT_M_TG1, m_TG1, em_TG1)

print('  Isolate TG-1:\n')
print('  T_opt =', T_optTG1, '+\-', Err_ToptTG1, '\n')


    # Isolate N-1

a_N1 = 0.0008722
T_L_N1 = 0
T_M_N1 = 35.586
m_N1 = 0.8964
ea_N1 = 0.0008244
eT_L_N1 = 0
eT_M_N1 = 2.6174
em_N1 = 0.2155

T_optN1 = ToptBriere(a_N1, T_L_N1, T_M_N1, m_N1)
Err_ToptN1 = Err_ToptBriere(a_N1, ea_N1, T_L_N1, eT_L_N1, T_M_N1, eT_M_N1, m_N1, em_N1)

print('  Isolate N-1:\n')
print('  T_opt =', T_optN1, '+\-', Err_ToptN1, '\n')


    # Isolate P.cor-441

a_PC44 = 0.002133
T_L_PC44 = 0
T_M_PC44 = 31.6897
m_PC44 = 1.1769
ea_PC44 = 0.001593
eT_L_PC44 = 0
eT_M_PC44 = 1.5273
em_PC44 = 0.3373

T_optPC44 = ToptBriere(a_PC44, T_L_PC44, T_M_PC44, m_PC44)
Err_ToptPC44 = Err_ToptBriere(a_PC44, ea_PC44, T_L_PC44, eT_L_PC44, T_M_PC44, eT_M_PC44, m_PC44, em_PC44)

print('  Isolate P.cor-441:\n')
print('  T_opt =', T_optPC44, '+\-', Err_ToptPC44, '\n')


    # Isolate P.cor-223

a_PC22 = 0.002133
T_L_PC22 = -9.03750
T_M_PC22 = 31.6945
m_PC22 = 1.4851637
ea_PC22 = 0.001593
eT_L_PC22 = 7
eT_M_PC22 = 2.11
em_PC22 = 0.71

T_optPC22 = ToptBriere(a_PC22, T_L_PC22, T_M_PC22, m_PC22)
Err_ToptPC22 = Err_ToptBriere(a_PC22, ea_PC22, T_L_PC22, eT_L_PC22, T_M_PC22, eT_M_PC22, m_PC22, em_PC22)

print('  Isolate P.cor-223:\n')
print('  T_opt =', T_optPC22, '+\-', Err_ToptPC22, '\n')
