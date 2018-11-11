# Sky Hoffert
# November 9, 2018

import cmath
from math import *

##########################################################################################################################
# functions
##########################################################################################################################
def reflection_coeff(Z_L, Z_0):
    return (Z_L - Z_0) / (Z_L + Z_0)

def parallel(Z_L, Z_0):
    return 1.0 / (1.0/Z_L + 1.0/Z_0)

##########################################################################################################################
# math constants not already defined
##########################################################################################################################
CONST_C = 3.0e8

##########################################################################################################################
# global values
##########################################################################################################################
f_low = 3.55e9
f_high = 3.7e9
f_center = (f_low + f_high)/2

##########################################################################################################################
# main program
##########################################################################################################################
print('================ Various ================')
print('Center Frequency: {:.3E}'.format(f_center))

##########################################################################################################################
# front end filter
##########################################################################################################################
print('\n================ Front End Filter ================')
# Using the three lambda/4 stub network
epsilon_r = 3
epsilon_eff = (epsilon_r + 1)/2 # for RO3003 material
lambda_0 = CONST_C / f_center
lambda_actual = lambda_0 / sqrt(epsilon_eff)

print('lambda_actual = {:.3E}'.format(lambda_actual))

h = 0.13 # mm
W = 0.32
W_prime = W
psi = (14 + 8/epsilon_r)/11 * (4* h / W_prime)
Z_c = 42.4 / sqrt(epsilon_r + 1) * log(1 + (4 * h / W_prime) * (psi + sqrt(psi**2 + (1 + 1/epsilon_r)/2 * pi**2)))

print('W = {:.3f}'.format(W))
print('Z_c = {:.3f}'.format(Z_c))

Z_stub = 1j * Z_c * tan(2 * pi * (lambda_actual/4)/lambda_actual)
print('Z_stub(f_c) = {:.3E}'.format(Z_stub))

f_off = 2.6e9
#f_off = f_center
Z_stub_2 = 1j * Z_c * tan(2 * pi * (lambda_actual/4)/(CONST_C / f_off / sqrt(epsilon_eff)))
print('Z_stub({:.3E} Hz) = {:.3E}'.format(f_off, Z_stub_2))

beta = 2 * pi / lambda_actual

# working right to left
Z_1 = parallel(Z_stub_2, Z_c)
Gamma_1L = (Z_1 - Z_c) / (Z_1 + Z_c)
Z_2 = Z_c * (1 + Gamma_1L * e**(-2 * 1j * beta * lambda_actual/4)) / (1 - Gamma_1L * e**(-2 * 1j * beta * lambda_actual/4))
Z_3 = parallel(Z_stub_2, Z_2)
Gamma_3L = (Z_3 - Z_c) / (Z_3 + Z_c)
Z_4 = Z_c * (1 + Gamma_3L * e**(-2 * 1j * beta * lambda_actual/4)) / (1 - Gamma_3L * e**(-2 * 1j * beta * lambda_actual/4))
Z_5 = parallel(Z_stub_2, Z_4)
Z_in = Z_5
print('Z_in = {:.3f}'.format(Z_in))

# s parameters for parallel impedance
s_21 = 2 * Z_c / (Z_in + 2 * Z_c)
s_12 = s_21
print('s_21 = {:.3f}'.format(s_21))
s_11 = (Z_in - Z_c) / (Z_in + Z_c)
#s_11 = Z_in / (Z_in + 2 * Z_c)
s_22 = s_11
print('s_11 = {:.3f}'.format(s_11))

Gamma_S = (Z_in - Z_c) / (Z_in + Z_c)
Gamma_L = (Z_c - Z_in) / (Z_c + Z_in)

print('Gamma_S = {:.3f}'.format(Gamma_S))
print('Gamma_L = {:.3f}'.format(Gamma_L))

TPG = ((1 - abs(Gamma_S)**2) * abs(s_21)**2 * (1 - abs(Gamma_L)**2)) / abs( (1 - s_11 * Gamma_S) * (1 - s_22 * Gamma_L) - s_12 * s_21 * Gamma_S * Gamma_L )**2
print('TPG ({:.3E} GHz) = {:.3f}'.format(f_off, TPG))

NF_fefilter = 10**(0.1/10)
G_fefilter = 1/NF_fefilter

##########################################################################################################################
# first amp
##########################################################################################################################
print('\n================ First Amp ================')

V_CC = 13.0 # V
I_C = 10.0e-3 # A
V_CE = 2.0 # V
V_BE = 0.7 # V
beta_min = 160.0
beta_max = 400.0
beta = 280.0 # 160 - 400

R_C = (V_CC - V_BE) / I_C
print('R_C = {:.3E}'.format(R_C))

R_C = 1100 # Ohm

R_1_smallest = beta_min * R_C
print('R_1_smallest = {:.3E}'.format(R_1_smallest))

R_1 = 20.5e3
I_C = (V_CC - V_BE) / (R_C + (R_1 + R_C)/beta_min)
I_B = ((V_CC - V_BE) - R_C * I_C) / (R_C + R_1)
V_CE = V_CC - (I_C + I_B) * R_C
P = V_CE * I_C

print('Final values:')
print('R_C = {:.3E}'.format(R_C))
print('R_1 = {:.3f}'.format(R_1))
print('I_B = {:.3E}'.format(I_B))
print('I_C = {:.3E}'.format(I_C))
print('V_CE = {:.3f} V'.format(V_CE))
print('Power = {:.3E} W'.format(P))

print('Blocking Cap:')
C_block = 100.0e-12
Z_block = -1j / (2 * pi * f_center * C_block)
print('C_block = {:.3E}'.format(C_block))
print('Z_block = {}'.format(Z_block))

# assume approximately unilateral (for now)
s_12 = 0.01

# assume V_A (early voltage value)
V_A = 100

s_21_2_dB = 20 # dB
s_21_2 = 10**(s_21_2_dB / 10)
s_21 = sqrt(s_21_2)
print('s_21^2 = {:.3f}'.format(s_21_2))
print('s_21 = {:.3f}'.format(s_21))

print('Calculating s_11:')
C_BE = 400.0e-15
r_i = beta / (40 * 0.04)
print('r_i = {:.3E}'.format(r_i))
Z_BE = -1j / (2 * pi * f_center * C_BE) + r_i
print('Z_BE = {} = {}'.format(Z_BE, cmath.polar(Z_BE)))
print('ang = {:.6f}'.format(degrees(cmath.polar(Z_BE)[1])))

print('Calculating s_22:')
r_o = (V_A + V_CE) / I_C
print('r_o = {:.3E}'.format(r_o))
C_CE = 268.0e-15
Z_CE = parallel(-1j / (2 * pi * f_center * C_CE), r_o)
print('Z_CE = {} = {}'.format(Z_CE, cmath.polar(Z_CE)))
print('ang = {:.6f}'.format(degrees(cmath.polar(Z_CE)[1])))

print('Calculating s_21:')
Z_CB = parallel(Z_BE, Z_CE) * s_21
print('Z_CB = {} = {}'.format(Z_CB, cmath.polar(Z_CB)))
print('ang = {:.6f}'.format(degrees(cmath.polar(Z_CB)[1])))

print('Normalize values to s_21:')
factor = (Z_CB / cmath.polar(Z_CB)[0])
s_11 = Z_BE / cmath.polar(Z_CB)[0]
s_22 = Z_CE / cmath.polar(Z_CB)[0]
s_21 *= factor
print('s_21 = {:.3f} = {:.3f} /_ {:.3f}'.format(s_21, abs(s_21), degrees(cmath.polar(s_21)[1])))
print('s_11 = {:.3f} = {:.3f} /_ {:.3f}'.format(s_11, abs(s_11), degrees(cmath.polar(s_11)[1])))
print('s_22 = {:.3f} = {:.3f} /_ {:.3f}'.format(s_22, abs(s_22), degrees(cmath.polar(s_22)[1])))

print('Stability')
Gamma_S = s_11.conjugate()
Gamma_L = s_22.conjugate()
Gamma_in = s_11 + (s_12 * s_21 * Gamma_L) / (1 - s_22 * Gamma_L)
Gamma_out = s_22 + (s_12 * s_21 * Gamma_S) / (1 - s_11 * Gamma_S)
delta = s_11 * s_22 - s_12 * s_21
K = (1 + abs(delta)**2 + abs(s_11)**2 + abs(s_22)**2) / (2 * abs(s_21 * s_12))
print('|delta| = {:.3f}'.format(abs(delta)))
print('K = {:.3f}'.format(K))
print('|Gamma_S| = {}'.format(abs(Gamma_S)))
print('|Gamma_L| = {}'.format(abs(Gamma_L)))
print('|Gamma_in| = {}'.format(abs(Gamma_in)))
print('|Gamma_out| = {}'.format(abs(Gamma_out)))

print('Gain:')
G_0 = abs(s_21)**2
G_1 = (1 - abs(Gamma_S)**2) / abs(1 - s_11 * Gamma_S)**2
G_2 = (1 - abs(Gamma_L)**2) / abs(1 - s_22 * Gamma_L)**2
print(G_0)
print(G_1)
print(G_2)
G_TU = G_1 * G_0 * G_2
print('G_TU = {} = {:.3f} dB'.format(G_TU, 10 * log10(G_TU)))

NF_amp1 = 10**(((0.7 + 0.47) / 2)/10)
G_amp1 = G_TU

##########################################################################################################################
# second amp
##########################################################################################################################
print('\n================ Second Amp ================')

NF_amp2 = NF_amp1
G_amp2 = G_TU

##########################################################################################################################
# Mixer/LO
##########################################################################################################################
print('\n================ Mixer ================')

NF_mixer = 10**(5.8/10)
G_mixer = 1/NF_mixer

##########################################################################################################################
# post-mixer filter
##########################################################################################################################
print('\n================ Post Mixer Filter ================')

f_down_center = 400.0e6
omega_center = 2 * pi * f_down_center
f_down_cutoff = 450.0e6
omega_cutoff = 2 * pi * f_down_cutoff
K = 11
Z_0 = 50

f = 500.0e6

G_f = (1 + (2 * pi * f / omega_cutoff)**(2 * K))**(-1)
G_center = (1 + (omega_center / omega_cutoff)**(2 * K))**(-1)

print('G_f ({:.3E} Hz) = {:.3f} = {:.3f} dB'.format(f, G_f, 10*log10(G_f)))
print('G_center {:.3E} Hz = {:.3f} = {:.3f} dB'.format(f_down_center, G_center, 10*log10(G_center)))

# filter values
for k in range(1, K+1):
    p_k = 2 * sin((2*k-1)*pi / (2 * K))
    p_k /= omega_cutoff
    if k % 2:
        p_k = p_k / Z_0
        print('p_{} = {:.3E} F'.format(k, p_k))
    else:
        p_k = p_k * Z_0
        print('p_{} = {:.3E} H'.format(k, p_k))

NF_pmfilter = 1/G_center
G_pmfilter = G_center

##########################################################################################################################
# Bias Tee
##########################################################################################################################
print('\n================ Bias Tee ================')

LC = 1 / omega_center**2

print('LC = {:.3E}'.format(LC))

C = 10.0e-12
L = LC / C

print('C = {:.3E}'.format(C))
print('L = {:.3E}'.format(L))

NF_biastee = 1
G_biastee = 1

##########################################################################################################################
# Tx Line to Base Station
##########################################################################################################################
print('\n================ TX line to base station ================')

Z_in = 75.0 # ohm
length = 10000 # m
attenuation = 10**(11.749/10) / 100.0 # per m @ 400 MHz

NF_cable = attenuation * length

##########################################################################################################################
# System
##########################################################################################################################
print('\n================ System ================')
print('Noise Factors:')
print('Front End Filter = {:.3f}'.format(NF_fefilter))
print('Amp 1 = {:.3f}'.format(NF_amp1))
print('Amp 2 = {:.3f}'.format(NF_amp2))
print('Mixer = {:.3f}'.format(NF_mixer))
print('Post Mix Filter = {:.3f}'.format(NF_pmfilter))
print('Bias tee = {:.3f}'.format(NF_biastee))
print('Cable = {:.3f}'.format(NF_cable))

# TODO - this calculation is wrong
NF_total = NF_fefilter + (NF_amp1-1)/(G_fefilter) + (NF_amp2-1)/(G_fefilter*G_amp1) + (NF_mixer-1)/(G_fefilter*G_amp1*G_amp2) + (NF_pmfilter-1)/(G_fefilter*G_amp1*G_amp2*G_mixer) + (NF_biastee-1)/(G_fefilter*G_amp1*G_amp2*G_mixer*G_pmfilter) + (NF_cable-1)/(G_fefilter*G_amp1*G_amp2*G_mixer*G_pmfilter*G_biastee)

print('Total NF = {:.3f} = {:.3f} dB'.format(NF_total, 10*log10(NF_total)))