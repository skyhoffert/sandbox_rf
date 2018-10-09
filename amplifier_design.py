# Sky Hoffert
# October 9, 2018

from cmath import *
from math import *
    
def display(name, val):
    print('{}: {:.4f}'.format(name, val))

# Given Problem Values
f = 2e9
Z_0 = 50
Z_in = 50
Z_out = 50
I_C = 10e-3
V_CE = 8
V_CC = 12
V_BE = 0.65
beta = 100

s_11 = rect(0.61, radians(165))
s_21 = rect(3.72, radians(57))
s_12 = rect(0.05, radians(42))
s_22 = rect(0.45, radians(-48))
'''
# Example values
s_11 = rect(0.30, radians(160))
s_21 = rect(6.10, radians(65))
s_12 = rect(0.03, radians(62))
s_22 = rect(0.40, radians(-38))
'''

print('==================== Stability ==========================')

delta = s_11 * s_22 - s_12 * s_21
K = (1 + pow((abs(delta)), 2) - pow(abs(s_11),2) - pow(abs(s_22),2)) / (2 * abs(s_12 * s_21))

display('delta', delta)
display('K', K)

print('==================== Gain ===============================')

MAG = abs(s_21 / s_12) * (K - pow(pow(K,2) - 1, 1/2))

display('MAG', MAG)
display('MAG (dB)', 10*log10(MAG))

print('==================== Matching ===========================')

B_1 = 1 - pow(abs(delta),2) + pow(abs(s_11),2) - pow(abs(s_22),2)
B_2 = 1 - pow(abs(delta),2) - pow(abs(s_11),2) + pow(abs(s_22),2)
M = s_11 - delta * s_22.conjugate()
N = s_22 - delta * s_11.conjugate()

display('B_1', B_1)
display('B_2', B_2)
display('M', M)
display('N', N)

Gamma_S = (B_1 - pow(pow(B_1,2) - 4 * pow(abs(M),2),1/2)) / (2 * M)
Gamma_L = (B_2 - pow(pow(B_2,2) - 4 * pow(abs(N),2),1/2)) / (2 * N)

print('Gamma_S: ', Gamma_S)
print('Gamma_L: ', Gamma_L)

Z_S = Z_0 * (1 + Gamma_S) / (1 - Gamma_S)
Z_L = Z_0 * (1 + Gamma_L) / (1 - Gamma_L)
display('Z_S', Z_S)
display('Z_L', Z_L)

print('==================== Confirming Gain ====================')

G_T = ((1 - pow(abs(Gamma_S),2)) * pow(abs(s_21),2) * (1 - pow(abs(Gamma_L),2))) / (pow(abs((1 - s_11*Gamma_S) * (1 - s_22*Gamma_L) - (s_12*s_21*Gamma_S*Gamma_L)),2))

display('G_T', G_T)
display('G_T (dB)', 10*log10(G_T))

print('==================== IMN ================================')

IMN_Y_0 = 1 / Z_0
IMN_Y_in = 1 / Z_in
IMN_G_in = IMN_Y_in.real
IMN_B_in = IMN_Y_in.imag
IMN_Z_L = Z_S.conjugate()

IMN_Gamma = (IMN_Z_L - Z_0) / (IMN_Z_L + Z_0)

display('IMN_Gamma', IMN_Gamma)

l = 0.001
increment = 0.001
wiggle = 0.0005
end = 1
while l < end:
    IMN_factor = IMN_Gamma * rect(1, -4*pi*l)
    IMN_Z_1 = Z_0 * (1 + IMN_factor) / (1 - IMN_factor)
    IMN_Y_1 = rect(1 / polar(IMN_Z_1)[0], -polar(IMN_Z_1)[1])
    IMN_G_1 = IMN_Y_1.real
    if abs(IMN_G_1 - IMN_G_in) < wiggle:
        break
    l += increment

display('l', l)
display('IMN_Y_1', IMN_Y_1)

l2_oc = 0.001
while l2_oc < end:
    IMN_Y_OC = IMN_Y_0 * tan(2 * pi * l2_oc) * 1j
    if abs(IMN_Y_OC.imag + IMN_Y_1.imag) < wiggle:
        break
    l2_oc += increment

display('l2_oc', l2_oc)
display('IMN_Y_OC', IMN_Y_OC)
    
l2_sc = 0.001
while l2_sc < end:
    IMN_Y_SC = IMN_Y_0 * 1/tan(2 * pi * l2_sc) * -1j
    if abs(IMN_Y_SC.imag + IMN_Y_1.imag) < wiggle:
        break
    l2_sc += increment

display('l2_sc', l2_sc)
display('IMN_Y_SC', IMN_Y_SC)

print('==================== OMN ================================')

OMN_Y_0 = 1 / Z_0
OMN_Y_in = 1 / Z_out
OMN_G_in = OMN_Y_in.real
OMN_B_in = OMN_Y_in.imag
OMN_Z_L = Z_L

OMN_Gamma = (OMN_Z_L - Z_0) / (OMN_Z_L + Z_0)

display('OMN_Gamma', OMN_Gamma)

l = 0.001
increment = 0.001
wiggle = 0.0005
end = 1
while l < end:
    OMN_factor = OMN_Gamma * rect(1, -4*pi*l)
    OMN_Z_1 = Z_0 * (1 + OMN_factor) / (1 - OMN_factor)
    OMN_Y_1 = rect(1 / polar(OMN_Z_1)[0], -polar(OMN_Z_1)[1])
    OMN_G_1 = OMN_Y_1.real
    if abs(OMN_G_1 - OMN_G_in) < wiggle:
        break
    l += increment

display('l', l)
display('OMN_Y_1', OMN_Y_1)

l2_oc = 0.001
while l2_oc < end:
    OMN_Y_OC = OMN_Y_0 * tan(2 * pi * l2_oc) * 1j
    if abs(OMN_Y_OC.imag + OMN_Y_1.imag) < wiggle:
        break
    l2_oc += increment

display('l2_oc', l2_oc)
display('OMN_Y_OC', OMN_Y_OC)
    
l2_sc = 0.001
while l2_sc < end:
    OMN_Y_SC = OMN_Y_0 * 1/tan(2 * pi * l2_sc) * -1j
    if abs(OMN_Y_SC.imag + OMN_Y_1.imag) < wiggle:
        break
    l2_sc += increment

display('l2_sc', l2_sc)
display('OMN_Y_SC', OMN_Y_SC)


print('==================== Biasing ============================')

R_C_R_E = (V_CC - V_CE) / (I_C)
R_C = R_C_R_E / 2
R_E = R_C

display('R_C', R_C)
display('R_E', R_E)

V_BB = I_C * R_E + V_BE

display('V_BB', V_BB)

R_1_R_2 = V_CC / (I_C / 10)

display('R_1 + R_2', R_1_R_2)

R_2 = V_BB * R_1_R_2 / V_CC
R_1 = R_1_R_2 - R_2
R_B = R_1 * R_2 / (R_1_R_2)

display('R_1', R_1)
display('R_2', R_2)

I_C_calc = (V_BB - V_BE) / ((R_B/beta) + R_E * beta / (beta + 1))

display('Actual I_C', I_C_calc)

while abs(I_C_calc - I_C) > 0.0001:
    R_E *= 0.99
    I_C_calc = (V_BB - V_BE) / ((R_B/beta) + R_E * beta / (beta + 1))

display('New I_C', I_C_calc)
display('New R_E', R_E)

R_C = R_C_R_E - R_E

display('New R_C', R_C)

V_CE_calc = V_CC - I_C * R_C - I_C * beta / (beta + 1) * R_E

display('Actual V_CE', V_CE_calc)


