# Sky Hoffert
# November 9, 2018

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
Z_stub_2 = 1j * Z_c * tan(2 * pi * (lambda_actual/4)/(CONST_C / f_off / sqrt(epsilon_eff)))
print('Z_stub({:.3E} Hz) = {:.3E}'.format(f_off, Z_stub_2))

# working right to left
Z_1 = Z_stub_2
Z_2 = Z_c + Z_1
Z_3 = parallel(Z_stub_2, Z_2)
Z_4 = Z_c + Z_3
Z_5 = parallel(Z_stub_2, Z_4)
Z_in = Z_5
print('Z_in = {:.3f}'.format(Z_in))

# working right to left
#Gamma_1 = reflection_coeff(parallel(50, Z_stub_2), 50)
#print('Gamma_1 = {:.3f}'.format(Gamma_1))
#Z_1 = Z_c * (1 + Gamma_1 * e**(-1j * 4 * pi * (lambda_actual/4)/(CONST_C / 2.6e9 / sqrt(epsilon_eff)))/(1 - Gamma_1 * e**(-1j * 4 * pi * (lambda_actual/4)/(CONST_C / 2.6e9 / #sqrt(epsilon_eff)))))
#print('Z_1 = {:.3f}'.format(Z_1))
#Gamma_2 = reflection_coeff(parallel(Z_1, Z_stub_2), 50)
#Z_2 = Z_c * (1 + Gamma_2 * e**(-1j * 4 * pi * (lambda_actual/4)/(CONST_C / 2.6e9 / sqrt(epsilon_eff)))/(1 - Gamma_2 * e**(-1j * 4 * pi * (lambda_actual/4)/(CONST_C / 2.6e9 / #sqrt(epsilon_eff)))))
##print('Z_2 = {:.3f}'.format(Z_2))
#Z_in = parallel(Z_2, Z_stub_2)
#print('Z_in = {:.3f}'.format(Z_in))
#Z_seen = 50 + 50**2 / Z_in

# s parameters for parallel impedance
s_21 = 2 * Z_c / (2 * Z_in + Z_c)
s_12 = s_21
print('s_21 = {:.3f}'.format(s_21))
#s_11 = (Z_in - Z_c) / (Z_in + Z_c)
s_11 = -Z_c / (2 * Z_in + Z_c)
s_22 = s_11
print('s_11 = {:.3f}'.format(s_11))

Gamma_S = (Z_in - Z_c) / (Z_in + Z_c)
Gamma_L = (Z_c - Z_in) / (Z_c + Z_in)

print('Gamma_S = {:.3f}'.format(Gamma_S))
print('Gamma_L = {:.3f}'.format(Gamma_L))

TPG = ((1 - abs(Gamma_S)**2) * abs(s_21)**2 * (1 - abs(Gamma_L)**2)) / abs( (1 - s_11 * Gamma_S) * (1 - s_22 * Gamma_L) - s_12 * s_21 * Gamma_S * Gamma_L )**2
print('TPG ({:.3E} GHz) = {:.3f}'.format(f_off, TPG))

NF_fefilter = 0.1 # dB

##########################################################################################################################
# first amp
##########################################################################################################################
print('\n================ First Amp ================')


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

NF_pmfilter = 10*log10(G_center) * -1

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

NF_biastee = 0

##########################################################################################################################
# Tx Line to Base Station
##########################################################################################################################
print('\n================ TX line to base station ================')

Z_in = 75.0 # ohm
attenuation = 11.749 / 100.0 # per m @ 400 MHz

##########################################################################################################################
# System
##########################################################################################################################
print('\n================ System ================')
print('Noise Factors:')
print('Front End Filter = {:.3f}'.format(NF_fefilter))
print('Post Mix Filter = {:.3f}'.format(NF_pmfilter))
print('Bias tee = {:.3f}'.format(NF_biastee))