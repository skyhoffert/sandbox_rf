# Sky Hoffert
# October 9, 2018

from cmath import *
from math import *
    
def display(name, val):
    print('{}: {:.4f}'.format(name, val))

# Given Problem Values
f = 1e9
Z_0 = 50
Z_in = 50
Z_out = 50
TPG_Desired_dB = 15
TPG_Desired = pow(10, TPG_Desired_dB/10)


s_11 = rect(0.73, radians(175))
s_21 = rect(4.45, radians(65))
s_12 = rect(0.05, radians(77))
s_22 = rect(0.21, radians(-80))
'''
# Example values
s_11 = rect(0.80, radians(120))
s_21 = rect(4.00, radians(60))
s_12 = rect(0.02, radians(62))
s_22 = rect(0.20, radians(-30))
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

print('==================== Unilateral Assumption ===========================')

G_1 = (1)/(1 - pow(abs(s_11),2))
G_0 = pow(abs(s_21),2)
G_2 = (1)/(1 - pow(abs(s_22),2))
G_TU_max = G_1 * G_0 * G_2

display('G_1', G_1)
display('G_1 (dB)', 10*log10(G_1))
display('G_0', G_0)
display('G_0 (dB)', 10*log10(G_0))
display('G_2', G_2)
display('G_2 (dB)', 10*log10(G_2))
display('G_TU_max', G_TU_max)
display('G_TU_max(dB)', 10*log10(G_TU_max))

u = abs(s_11 * s_12 * s_21 * s_22) / abs((1 - pow(abs(s_11),2)) * (1 - pow(abs(s_22),2)))

display('u (Unilateral Figure of Merit)', u)

err_low = 1/pow(abs(1 + u), 2)
err_high = 1/pow(abs(1 - u), 2)

display('err_low', err_low)
display('err_low (dB)', 10*log10(err_low))
display('err_high', err_high)
display('err_high (dB)', 10*log10(err_high))

# adjust G_1 to meet desired TPG
diff = (10*log10(G_TU_max) - TPG_Desired_dB)
display('diff', diff)
G_1 = pow(10, (10*log10(G_1) - diff)/10)
display('New G_1', G_1)
display('New G_1 (dB)', 10*log10(G_1))

g_1 = G_1 * (1 - pow(abs(s_11),2))

display('g_1', g_1)

C_1 = (g_1 * s_11.conjugate()) / (1 - pow(abs(s_11),2) * (1 - g_1))
r_1 = (pow(1 - g_1, 1/2) * (1 - pow(abs(s_11),2))) / (1 - pow(abs(s_11),2) * (1 - g_1))

display('C_1', C_1)
print('C_1 (polar): ', polar(C_1))
display('r_1', r_1)

Gamma_S = (abs(C_1) - r_1) * C_1 / abs(C_1)
Gamma_L = s_22.conjugate()

display('Gamma_S', Gamma_S)
print('Gamma_S (polar): ', polar(Gamma_S))
display('Gamma_L', Gamma_L)

print('==================== Matching ===========================')

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
OMN_Z_L = Z_L.conjugate()

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


