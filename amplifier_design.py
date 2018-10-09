# Sky Hoffert
# October 9, 2018

from cmath import *
from math import *
    
def display(name, val):
    print('{}: {:.3f}'.format(name, val))

# Given Problem Values
f = 2e9
Z_0 = 50
Z_in = 50
Z_out = 50
I_C = 10e-3
V_CE = 8
V_CC = 12

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

delta = s_11 * s_22 - s_12 * s_21

K = (1 + pow((abs(delta)), 2) - pow(abs(s_11),2) - pow(abs(s_22),2)) / (2 * abs(s_12 * s_21))

MAG = abs(s_21 / s_12) * (K - pow(pow(K,2) - 1, 1/2))

B_1 = 1 - pow(abs(delta),2) + pow(abs(s_11),2) - pow(abs(s_22),2)
B_2 = 1 - pow(abs(delta),2) - pow(abs(s_11),2) + pow(abs(s_22),2)
M = s_11 - delta * s_22.conjugate()
N = s_22 - delta * s_11.conjugate()

Gamma_S = (B_1 - pow(pow(B_1,2) - 4 * pow(abs(M),2),1/2)) / (2 * M)
Gamma_L = (B_2 - pow(pow(B_2,2) - 4 * pow(abs(N),2),1/2)) / (2 * N)

Z_S = Z_0 * (1 + Gamma_S) / (1 - Gamma_S)
Z_L = Z_0 * (1 + Gamma_L) / (1 - Gamma_L)

# display values
print('======== Stability ========')
display('delta', delta)
display('K', K)
print('======== Gain ========')
display('MAG', MAG)
display('MAG (dB)', 10*log10(MAG))
print('======== Matching ========')
display('B_1', B_1)
display('B_2', B_2)
display('M', M)
display('N', N)
print('Gamma_S: ', Gamma_S)
print('Gamma_L: ', Gamma_L)
display('Z_S', Z_S)
display('Z_L', Z_L)

