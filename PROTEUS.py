#==============================================================================================================
#
#                Photochemical and RadiatiOn Transport model for Extensive USe (PROTEUS)
#
#==============================================================================================================
#
# 2 July 2020, Yuki Nakamura
#
#    * Users are needed to install 'tkinter' to run this GUI.
#      Python 3 is recommended.
#
#    < for Python 3 >
#      'tkinter' is treated as 'tkinter' in this version. So you need to import 'tkinter as tk'.
#

import tkinter as tk #GUI library
from tkinter import *
import re #Regular expression
import os #operating system interface
import matplotlib.pyplot as plt #plot library
import numpy as np #mathematical library
import webbrowser #To open weblink

version = 1.0

# define reaction_rate_list
reaction_rate_list = []
doi_list = []
Planet_list = []
iplnt = 0

# reaction & rate list ####################################################################################################
# NOTES :
#   * " reaction : rate coefficient @ reference # label "
#   * Reaction and rate should be separated by ":" , rate and reference/label should be separeted by @/#, respectively.
#   * In reaction expression, each species and array should be separeted by at least 1 space.
#   * Expression of chemical species
#       Normal rule   : e.g.) Carbon dioxide = CO2, water cluster = (H2O)2
#       Excited states: **** to be written ****
#         
#       Electron      : e-
#       Positive ions : e.g.) CO2+, O(2D)+, H+(H2O)2, CO2++
#       Negative ions : e.g.) CO3-, NO3-(H2O)3
#
#   * Temperature expression
#       Neutral temperature : Tn
#       Ion temperature     : Ti
#       Electron temperature: Te
#
#   * Altitude dependent: altitude should be expressed as 'h' in [km]
#

#--------------------------------------------------------------------------------------------------------------------------
#
#                                                       Venus
#
#--------------------------------------------------------------------------------------------------------------------------
Planet_list.append(['Venus',len(reaction_rate_list)])
# reaction list
# neutral
# Krasnopolsky 2007
reaction_rate_list.append(" SO3 + H2O + H2O -> H2SO4 + H2O : 2.3e-43 * Tn * exp(6540/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" H2SO4 + H2O -> SO3 + H2O + H2O : 7.0e-14 * exp(-5170/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" SO3 + CO -> CO2 + SO2 : 1.0e-11 * exp(-13000/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" SO3 + OCS -> CO2 + SO2 + (SO)2 : 1.0e-11 * exp(-10000/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" (SO)2 + OCS -> CO + SO2 + S2 : 1.0e-20 # Krasnopolsky [2007] ")
reaction_rate_list.append(" SO + SO -> SO2 + S : 1.0e-12 * exp(-1700/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + SO2 -> SO + SO : 2.3e-11 * exp(-5200/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" CO + SO2 -> CO2 + SO : 4.5e-12 * exp(-24300/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" SO + CO2 -> SO2 + CO : 1.5e-11 * exp(-22000/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S3 -> S2 + S : 0.00023 * exp(0.15*h - 0.00119*h^2) # Krasnopolsky [2007] ") # photochemical reaction
reaction_rate_list.append(" S + S + M -> S2 + M : 1.0e-30 * (300/Tn)^2 # Krasnopolsky [2007] ")
reaction_rate_list.append(" S2 + M -> S + S + M : 2.2e-7 * exp(-50000/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + S2 + M -> S3 + M : 1.0e-30 * (300/Tn)^2 # Krasnopolsky [2007] ")
reaction_rate_list.append(" S3 + M -> S + S2 + M : 1.3e-6 * exp(-29800/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + CO + M -> OCS + M : 3.0e-33 * exp(-1000/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" OCS + M -> CO + S + M : 2.2e-7 * exp(-37300/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + OCS -> CO + S2 : 1.7e-11 * exp(-2800/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S2 + CO -> OCS + S : 1.0e-12 * exp(-17460/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + S3 -> S2 + S2 : 1.7e-10 * exp(-2800/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S2 + S2 -> S + S3 : 2.8e-11 * exp(-23000/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" CO + S3 -> OCS + S2 : 1.0e-11 * exp(-20000/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S2 + OCS -> CO + S3 : 2.5e-11 * exp(-25500/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + NO + M -> SNO + M : 3.0e-32 * exp(940/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + SNO -> S2 + NO : 5.0e-11 # Krasnopolsky [2007] ")
reaction_rate_list.append(" S2 + SNO -> S3 + NO : 1.0e-17 # Krasnopolsky [2007] ")
reaction_rate_list.append(" SH + SH -> H2S + S : 1.5e-11 # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + H2S -> SH + SH : 1.7e-10 * exp(-3620/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + SH -> S2 + H : 4.5e-11 # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + S2 -> SH + S : 5.3e-10 * exp(-8830/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + SH -> H2 + S : 2.5e-11 # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + H2 -> SH + H : 1.0e-10 * exp(-10080/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + OCS -> CO + SH : 1.2e-11 * exp(-1950/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" CO + SH -> OCS + H : 6.3e-14 * exp(-7780/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + HCl -> H2 + Cl : 1.7e-11 * exp(-1770/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" Cl + H2 -> HCl + H : 3.0e-11 * exp(-2270/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" H2S + Cl -> HCl + SH : 3.7e-11 * exp(210/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" SH + HCl -> H2S + Cl : 7.6e-12 * exp(-5750/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" Cl + SH -> HCl + S : 8.0e-11 * exp(210/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + HCl -> SH + Cl : 1.9e-10 * exp(-9380/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + SNO -> NO + SH : 4.0e-10 * exp(-340/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + SH + M -> H2S + M : 1.0e-30 * (300/Tn)^2 # Krasnopolsky [2007] ")
reaction_rate_list.append(" H2S + M -> SH + H + M : 1.9e-7 * exp(-44750/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + S3 -> SH + S2 : 1.2e-10 * exp(-1950/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" Cl + SO2 + M -> ClSO2 + M : 1.3e-34 * exp(940/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" ClSO2 + M -> Cl + SO2 + M : 7.0e-16 * exp(-10540/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" ClSO2 + ClSO2 -> SO2Cl2 + SO2 : 1.0e-12 # Krasnopolsky [2007] ")
reaction_rate_list.append(" SO2Cl2 + SO2 -> ClSO2 + ClSO2 : 1.0e-12 * exp(-11000/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" ClSO2 + Cl -> SO2 + Cl2 : 1.0e-12 # Krasnopolsky [2007] ")
reaction_rate_list.append(" ClSO2 + S -> SO2 + SCl : 1.0e-12 # Krasnopolsky [2007] ")
reaction_rate_list.append(" ClSO2 + H -> SO2 + HCl : 1.0e-11 # Krasnopolsky [2007] ")
reaction_rate_list.append(" SO2Cl2 + Cl -> ClSO2 + Cl2 : 1.0e-12 # Krasnopolsky [2007] ")
reaction_rate_list.append(" SO2Cl2 + S -> ClSO2 + SCl : 1.0e-12 # Krasnopolsky [2007] ")
reaction_rate_list.append(" SO2Cl2 + H -> ClSO2 + HCl : 1.0e-11 # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + Cl2 -> HCl + Cl : 8.0e-11 * exp(-416/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + Cl2 -> SCl + Cl : 2.8e-11 * exp(-300/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" Cl + SCl -> S + Cl2 : 2.8e-11 * exp(-650/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + SCl -> S2 + Cl : 1.0e-12 # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + SCl -> HCl + S : 1.0e-11 # Krasnopolsky [2007] ")
reaction_rate_list.append(" SH + Cl2 -> HSCl + Cl : 1.4e-11 * exp(-690/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" HSCl + SH -> H2S + SCl : 3.0e-12 * exp(-500/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" HSCl + S -> SH + SCl : 1.7e-13 # Krasnopolsky [2007] ")
reaction_rate_list.append(" HSCl + H -> H2 + SCl : 1.2e-13 * exp(-2770/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" HSCl + Cl -> HCl + SCl : 2.5e-13 * exp(-130/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" SH + SCl -> S2 + HCl : 6.0e-13 * exp(230/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" SH + OH -> H2O + S : 2.5e-12 # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + H2O -> OH + SH : 4.7e-11 * exp(-17700/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" OH + H2 -> H2O + H : 2.8e-12 * exp(-1800/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + H2O -> OH + H2 : 1.3e-11 * exp(-9420/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" OH + HCl -> H2O + Cl : 2.6e-12 * exp(-350/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" Cl + H2O -> OH + HCl : 2.0e-11 * exp(-8470/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" OH + H2S -> H2O + SH : 6.1e-12 * exp(-75/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" SH + H2O -> H2S + OH : 1.0e-11 * exp(-14160/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + H2S -> SH + H2 : 8.2e-11 * exp(-1470/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" SH + H2 -> H2S + H : 3.0e-11 * exp(-7930/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" CO + OH -> CO2 + H : 2.4e-13 # Krasnopolsky [2007] ")
reaction_rate_list.append(" H + CO2 -> CO + OH : 1.0e-10 * exp(-12400/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" OH + OCS -> CO2 + SH : 1.1e-13 * exp(-1200/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" SH + CO2 -> OH + OCS : 2.6e-13 * exp(-19360/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" SO + OH -> SO2 + H : 2.7e-11 * exp(335/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" SO2 + H -> SO + OH : 3.7e-9 * exp(-14350/Tn) # Krasnopolsky [2007] ")
reaction_rate_list.append(" S + OH -> SO + H : 6.6e-11 # Krasnopolsky [2007] ")
reaction_rate_list.append(" SO + H -> S + OH : 4.0e-10 * exp(-11200/Tn) # Krasnopolsky [2007] ")


#--------------------------------------------------------------------------------------------------------------------------
#
#                                                       Earth
#
#--------------------------------------------------------------------------------------------------------------------------
Planet_list.append(['Earth',len(reaction_rate_list)])
# reaction list

# photoionization
reaction_rate_list.append(" O2  + hv -> O2+    + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" O2  + hv -> O+(4S) + e- + O  : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" O   + hv -> O+(4S) + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" O   + hv -> O+(2D) + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" O   + hv -> O+(2P) + e-      : Photoionization @ Schunk and Nagy [2001] ")
#reaction_rate_list.append(" O   + hv -> O+(4P) + e-      : Photoionization @ Schunk and Nagy [2001] ") # NO LOSS for O+(4P)
reaction_rate_list.append(" N2  + hv -> N2+    + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" N2  + hv -> N+     + e- + N  : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" H   + hv -> H+     + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" H2  + hv -> H2+    + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" H2  + hv -> H+     + e- + H  : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" He  + hv -> He+    + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" Ar  + hv -> Ar+    + e-      : Photoionization @ Schunk and Nagy [2001] ")

## photodissociation
# Chaffin et al. [2017]
#reaction_rate_list.append(" H2O  + hv -> H      + OH      : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H2O  + hv -> H2     + O(1D)   : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O3   + hv -> O2     + O       : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O3   + hv -> O2     + O(1D)   : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O2   + hv -> O      + O       : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O2   + hv -> O      + O(1D)   : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H2   + hv -> H      + H       : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" OH   + hv -> O      + H       : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" OH   + hv -> O(1D)  + H       : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" HO2  + hv -> OH     + O       : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H2O2 + hv -> OH     + OH      : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H2O2 + hv -> HO2    + H       : Photodissociation # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H2O2 + hv -> H2O    + O(1D)   : Photodissociation # Chaffin et al. [2017] ")

# photodissociation (Tian et al. [2011])
reaction_rate_list.append(" O2       + hv -> O         + O(1D)       : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" O2       + hv -> O         + O           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" H2O      + hv -> H         + OH          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" O3       + hv -> O2        + O(1D)       : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" O3       + hv -> O2        + O           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" H2O2     + hv -> OH        + OH          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CO2      + hv -> CO        + O           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CO2      + hv -> CO        + O(1D)       : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" H2CO     + hv -> H2        + CO          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" H2CO     + hv -> HCO       + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" HCO      + hv -> H         + CO          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" HO2      + hv -> OH        + O           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH4      + hv -> ^1CH2     + H2          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH4      + hv -> ^3CH2     + H      + H  : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH4      + hv -> CH3       + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C2H6     + hv -> ^3CH2     + ^3CH2  + H2 : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C2H6     + hv -> CH4       + ^1CH2       : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C2H6     + hv -> C2H2      + H2     + H2 : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C2H6     + hv -> C2H4      + H      + H  : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C2H6     + hv -> C2H4      + H2          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C2H6     + hv -> CH3       + CH3         : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" HNO2     + hv -> NO        + OH          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" HNO3     + hv -> NO2       + OH          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" NO       + hv -> N(4S)     + O           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" NO2      + hv -> NO        + O           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + hv -> ^1CH2     + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" HNO      + hv -> NO        + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" NH3      + hv -> NH2       + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" N2H4     + hv -> N2H3      + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" NH       + hv -> N(4S)     + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" NH2      + hv -> NH        + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C2H2     + hv -> C2H       + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C2H2     + hv -> C2        + H2          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C2H4     + hv -> C2H2      + H2          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C2H4     + hv -> C2H2      + H      + H  : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C3H8     + hv -> C3H6      + H2          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C3H8     + hv -> C2H6      + ^1CH2       : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C3H8     + hv -> C2H4      + CH4         : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C3H8     + hv -> C2H5      + CH3         : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C3H6     + hv -> C2H2      + CH3    + H  : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C3H6     + hv -> CH2CCH2   + H2          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C3H6     + hv -> C2H4      + ^3CH2       : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C3H6     + hv -> C2H       + CH4    + H  : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + hv -> C         + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH2CO    + hv -> ^3CH2     + CO          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH3CHO   + hv -> CH3       + HCO         : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH3CHO   + hv -> CH4       + CO          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5CHO  + hv -> C2H5      + HCO         : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" C3H3     + hv -> C3H2      + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH3C2H   + hv -> C3H3      + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH3C2H   + hv -> C3H2      + H2          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH3C2H   + hv -> CH3       + C2H         : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH2CCH2  + hv -> C3H3      + H           : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH2CCH2  + hv -> C3H2      + H2          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" CH2CCH2  + hv -> C2H2      + ^3CH2       : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" HCN      + hv -> H         + CN          : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" HNCO     + hv -> H         + NCO         : Photodissociation # Tian et al. [2011] ")
reaction_rate_list.append(" NCO      + hv -> N         + CO          : Photodissociation # Tian et al. [2011] ")

## neutral

# Chaffin et al. [2017]
#reaction_rate_list.append(" O     + O   + M   -> O2   + M        : 5.4e-33 * (300/Tn)^3.25 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O     + O2  + N2  -> O3   + N2       : 5.0e-35 * exp(724/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O     + O2  + CO2 -> O3   + CO2      : 1.5e-33 * (300/Tn)^2.4 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O     + O3        -> O2   + O2       : 8.0e-12 * exp(-2060/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O     + CO  + M   -> CO2  + M        : 2.2e-33 * exp(-1780/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O(1D) + O2        -> O    + O2       : 3.2e-11 * exp(70/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O(1D) + O3        -> O2   + O2       : 1.2e-10 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O(1D) + O3        -> O    + O  + O2  : 1.2e-10 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O(1D) + H2        -> H    + OH       : 1.2e-10 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O(1D) + CO2       -> O    + CO2      : 7.5e-11 * exp(115/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O(1D) + H2O       -> OH   + OH       : 1.63e-10 * exp(60/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H2    + O         -> OH   + H        : 6.34e-12 * exp(-4000/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" OH    + H2        -> H2O  + H        : 9.01e-13 * exp(-1526/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H     + H   + CO2 -> H2   + CO2      : 1.6e-32 * (298/Tn)^2.27 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H     + OH  + CO2 -> H2O  + CO2      : 1.292e-30 * (300/Tn)^2 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H     + HO2       -> OH   + OH       : 7.2e-11 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H     + HO2       -> H2O  + O(1D)    : 1.6e-12 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H     + HO2       -> H2   + O2       : 3.45e-12 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H     + H2O2      -> HO2  + H2       : 2.8e-12 * exp(-1890/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H     + H2O2      -> H2O  + OH       : 1.7e-11 * exp(-1800/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H     + O2  + M   -> HO2  + M        : k0 = 8.8e-32 * (300/Tn)^1.3 && \
#                                                                   kinf = 7.5e-11 * (300/Tn)^(-0.2) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" H     + O3        -> OH   + O2       : 1.4e-10 * exp(-470/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O     + OH        -> O2   + H        : 1.8e-11 * exp(180/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O     + HO2       -> OH   + O2       : 3.0e-11 * exp(200/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" O     + H2O2      -> OH   + HO2      : 1.4e-12 * exp(-2000/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" OH    + OH        -> H2O  + O        : 1.8e-12 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" OH    + OH  + M   -> H2O2 + M        : k0 = 8.97e-31 * (300/Tn)^1 && \
#                                                                   kinf = 2.6e-11 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" OH    + O3        -> HO2  + O2       : 1.7e-12 * exp(-940/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" OH    + HO2       -> H2O  + O2       : 4.8e-11 * exp(250/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" OH    + H2O2      -> H2O  + HO2      : 1.8e-12 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" HO2   + O3        -> OH   + O2 + O2  : 1.0e-14 * exp(-490/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" HO2   + HO2       -> H2O2 + O2       : 3.0e-13 * exp(460/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" HO2   + HO2 + M   -> H2O2 + O2 + M   : 4.2e-33 * exp(920/Tn) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" CO    + OH  + M   -> CO2  + H  + M   : k0M = 1.5e-13 * (300/Tn)^(-0.6) && \
#                                                                   kinfM = 2.1e9 * (300/Tn)^(-6.1) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" CO    + OH  + M   -> HOCO + M        : k0 = 5.9e-33 * (300/Tn)^1.4 && \
#                                                                   kinf = 1.1e-12 * (300/Tn)^(-1.3) # Chaffin et al. [2017] ")
#reaction_rate_list.append(" HOCO  + O2        -> HO2  + CO2      : 2.0e-12 # Chaffin et al. [2017] ")
#reaction_rate_list.append(" CO2+  + H2        -> CO2  + H  + H   : 8.7e-10 # Chaffin et al. [2017] ") # it's only used to reproduce Chaffin+2017, otherwise remove it!
#reaction_rate_list.append(" CO    + OH        -> CO2  + H        : 1.5e-13 # Chaffin et al. [2017] ")

# Tian et al. [2011]
reaction_rate_list.append(" H2O      + O(1D)          -> OH       + OH               : 1.63e-10 * exp(60/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H2       + O(1D)          -> OH       + H                : 1.00e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" H2       + O              -> OH       + H                : 8.5e-20 * (Tn)^2.67 * exp(-3163/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H2       + OH             -> H2O      + H                : 9.01e-13 * (Tn/298)^2.41 * exp(-1527/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H        + O3             -> OH       + O2               : 1.4e-10 * exp(-470/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H        + O2      + M    -> HO2      + M                : k0 = 4.4e-32 * (300/Tn)^1.3 && \
                                                                                       kinf = 4.7e-11 * (300/Tn)^(0.2) # Tian et al. [2011] ")
reaction_rate_list.append(" H        + HO2            -> H2       + O2               : 0.08 * 8.1e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" H        + HO2            -> H2O      + O                : 0.02 * 8.1e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" H        + HO2            -> OH       + OH               : 0.9 * 8.1e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" OH       + O              -> H        + O2               : 2.2e-11 * exp(120/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" OH       + HO2            -> H2O      + O2               : 4.8e-11 * exp(250/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" OH       + O3             -> HO2      + O2               : 1.7e-12 * exp(-940/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HO2      + O              -> OH       + O2               : 3.0e-11 * exp(2000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HO2      + O3             -> OH       + O2     + O2      : 1.0e-14 * exp(-490/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HO2      + HO2            -> H2O2     + O2               : 3.5e-13 * exp(430/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HO2      + HO2     + M    -> H2O2     + O2      + M      : 1.7e-33 * exp(1000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H2O2     + OH             -> HO2      + H2O              : 1.8e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" O        + O       + M    -> O2       + M                : 5.2e-35 * exp(900/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" O        + O2      + M    -> O3       + M                : k0 = 6.0e-34 * (300/Tn)^2.4 && \
                                                                                       kinf = 1.0e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" O        + O3             -> O2       + O2               : 8.0e-12 * exp(-2060/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" OH       + OH             -> H2O      + O                : 2.3e-20 * (Tn)^(2.6) * exp(946/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" O(1D)    + M              -> O        + M                : 7.5e-11 * exp(115/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" O(1D)    + O2             -> O        + O2               : 3.12e-11 * exp(70/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CO       + OH      + M    -> CO2      + H       + M      : k0M = 1.5e-13 * (300/Tn)^(-0.6) && \
                                                                                       kinfM = 2.1e9 * (300/Tn)^(-6.1) # Chaffin et al. [2017] ")
reaction_rate_list.append(" CO       + O       + M    -> CO2      + M                : k0 = 1.7e-33 * exp(-1510/Tn) && \
                                                                                       kinf = 2.66e-14 * exp(-1459/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H        + CO      + M    -> HCO      + M                : k0 = 1.4e-34 * exp(-100/Tn) && \
                                                                                       kinf = 1.96e-13 * exp(-1366/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H        + HCO            -> H2       + CO               : 1.8e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" HCO      + HCO            -> H2CO     + CO               : 4.5e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" OH       + HCO            -> H2O      + CO               : 1.7e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" O        + HCO            -> H        + CO2              : 5.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" O        + HCO            -> OH       + CO               : 5.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" H2CO     + H              -> H2       + HCO              : 2.1e-16 * (Tn)^(1.62) * exp(-1090/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H        + H       + M    -> H2       + M                : 2.7e-31 * (Tn)^(-0.6) # Tian et al. [2011] ")
reaction_rate_list.append(" HCO      + O2             -> HO2      + CO               : 5.2e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" H2CO     + OH             -> H2O      + HCO              : 8.2e-12 * exp(40/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H        + OH      + M    -> H2O      + M                : 6.1e-26 * (Tn)^(-2.0) # Tian et al. [2011] ")
reaction_rate_list.append(" OH       + OH      + M    -> H2O2     + M                : k0 = 6.9e-31 * (Tn)^(-1.0) && \
                                                                                       kinf = 2.6e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" H2O2     + O              -> OH       + HO2              : 1.4e-12 * exp(-2000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH4      + OH             -> CH3      + H2O              : 4.16e-13 * (Tn/298)^(2.18) * exp(-1232/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH4      + O(1D)          -> CH3      + OH               : 1.13e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CH4      + O(1D)          -> H2CO     + H2               : 7.5e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" ^1CH2    + CH4            -> CH3      + CH3              : 6.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" ^1CH2    + O2             -> HCO      + OH               : 3.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" ^1CH2    + M              -> ^3CH2    + M                : 8.8e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + H2             -> CH3      + H                : 5.0e-14 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + CH4            -> CH3      + CH3              : 7.1e-12 * exp(-5051/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + O2             -> HCO      + OH               : 1.5e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + O2             -> H2CO     + OH               : 5.5e-13 * exp(-4500/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + OH             -> H2CO     + H2               : 3.76e-14 * (Tn)^(-0.12) * exp(209/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + O              -> H2CO     + H                : 1.10e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + O3             -> H2CO     + HO2              : 5.4e-12 * exp(-220/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + CH3     + M    -> C2H6     + M                : k0 = 3.5e-7 * (Tn)^(-7) * exp(-1390/Tn) && \
                                                                                       kinf = 1.5e-7 * (Tn)^(-1.18) * exp(-329/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + H       + M    -> CH4      + M                : k0 = 1.7e-24 * (Tn)^(-1.8) && \
                                                                                       kinf = 3.5e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + HCO            -> CH4      + CO               : 8.20e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + HNO            -> CH4      + NO               : 3.00e-14 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + H2CO           -> CH4      + HCO              : 1.3e-31 * (Tn)^(6.1) * exp(-990/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H        + NO      + M    -> HNO      + M                : 2.1e-32 * exp(300/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" N(4S)    + N(4S)   + M    -> N2       + M                : 8.3e-34 * exp(500/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" N(4S)    + O2             -> NO       + O                : 1.5e-11 * exp(-3600/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" N(4S)    + O3             -> NO       + O2               : 2.0e-16 # Tian et al. [2011] ")
reaction_rate_list.append(" N(4S)    + OH             -> NO       + H                : 6.08e-11 * (Tn/298)^(-0.69) * exp(-48.1/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" N(4S)    + NO             -> N2       + O                : 2.1e-11 * exp(100/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" NO       + O3             -> NO2      + O2               : 3.0e-12 * exp(-1500/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" NO       + O       + M    -> NO2      + M                : k0 = 9.0e-32 * (Tn/300)^(-1.5) && \
                                                                                       kinf = 3.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" NO       + HO2            -> NO2      + OH               : 3.5e-12 * exp(250/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" NO       + OH      + M    -> HNO2     + M                : k0 = 7.0e-31 * (Tn/300)^(-2.6) && \
                                                                                       kinf = 3.6e-11 * (Tn/300)^(-0.1) # Tian et al. [2011] ")
reaction_rate_list.append(" NO2      + O              -> NO       + O2               : 5.1e-12 * exp(210/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" NO2      + OH      + M    -> HNO3     + M                : k0 = 9.1e-32 * (Tn/300)^(-3.0) && \
                                                                                       kinf = 4.2e-11 * (Tn/300)^(-0.5) # Tian et al. [2011] ")
reaction_rate_list.append(" NO2      + H              -> NO       + OH               : 4.0e-10 * exp(-340/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HNO3     + OH             -> H2O      + NO2      + O     : 2.4e-14 * exp(460/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HCO      + NO             -> HNO      + CO               : 1.3e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" H        + HNO            -> H2       + NO               : 3.0e-11 * exp(-500.6/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" O        + HNO            -> OH       + NO               : 3.8e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" OH       + HNO            -> H2O      + NO               : 6.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" HNO2     + OH             -> H2O      + NO2              : 1.8e-11 * exp(-390/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH4      + O              -> CH3      + OH               : 1.15e-15 * (Tn)^(1.56) * exp(-4270/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" ^1CH2    + H2             -> CH3      + H                : 1.2e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" ^1CH2    + CO2            -> H2CO     + CO               : 1.0e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + O              -> HCO      + H                : 1.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + CO2            -> H2CO     + CO               : 3.9e-14 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H6     + OH             -> C2H5     + H2O              : 8.7e-12 * exp(-1070/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H6     + O              -> C2H5     + OH               : 1.66e-15 * (Tn)^(1.5) * exp(-2920/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H6     + O(1D)          -> C2H5     + OH               : 1.5e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + H              -> CH3      + CH3              : 1.25e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + O              -> CH3      + HCO      + H     : 1.70e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + OH             -> CH3      + HCO      + H2    : 1.10e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + HCO            -> C2H6     + CO               : 2.00e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + HNO            -> C2H6     + NO               : 3.00e-14 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + O2      + M    -> CH3 + CH3 + HCO + OH + M    : k0 = 1.5e-28 * (Tn/300)^(-3.0) && \
                                                                                       kinf = 8.0e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" NH3      + OH             -> NH2      + H2O              : 1.7e-12 * exp(-710/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" NH3      + O(1D)          -> NH2      + OH               : 2.5e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" NH2      + H       + M    -> NH3      + M                : 3.01e-30 # Tian et al. [2011] ")
reaction_rate_list.append(" NH2      + NO             -> N2       + H2O              : 3.8e-12 * exp(450/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" NH2      + NH2     + M    -> N2H4     + M                : 1.0e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" NH2      + O              -> NH       + OH               : 5.0e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" NH2      + O              -> HNO      + H                : 5.0e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" NH       + NO             -> N2       + O    + H         : 4.9e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" NH       + O              -> N(4S)    + OH               : 1.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" N2H4     + H              -> N2H3     + H2               : 9.9e-12 * exp(-1200/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" N2H3     + H              -> NH2      + NH2              : 2.7e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" N2H3     + N2H3           -> N2H4     + N2   + H2        : 6.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" NH       + H       + M    -> NH2      + M                : 6.0e-30 # Tian et al. [2011] ")
reaction_rate_list.append(" NH2      + HCO            -> NH3      + CO               : 1.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" NH       + HCO            -> NH2      + CO               : 1.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" ^1CH2    + O2             -> H2CO     + O                : 3.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + CH3            -> C2H4     + H                : 7.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + CH3            -> C3H8                        : 6.4e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C3H8     + OH             -> C3H7     + H2O              : 8.7e-12  * exp(-615/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C3H8     + O              -> C3H7     + OH               : 2.2e-11  * exp(-2250/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C3H8     + O(1D)          -> C3H7     + OH               : 1.3e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C3H7     + H              -> CH3      + C2H5             : 6.0e-13 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + ^3CH2          -> C2H2     + H   + H          : 1.8e-10 * exp(-400/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H2     + OH             -> CO       + CH3              : 1.91e-12 * exp(-233/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H2     + H       + M    -> C2H3     + M                : k0 = 3.3e-30 * exp(-740/Tn) && \
                                                                                       kinf = 1.4e-11 * exp(-1300/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H3     + H              -> C2H2     + H2               : 6.86e-11 * exp(-23/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H3     + H2             -> C2H4     + H                : 1.57e-20 * (Tn)^(2.56) * exp(-2529/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H3     + CH4            -> C2H4     + CH3              : 2.4e-24 * (Tn)^(4.02) * exp(-2754/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H3     + C2H6           -> C2H4     + C2H5             : 1.0e-21 * (Tn)^(3.3) * exp(-5285/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H4     + OH             -> H2CO     + CH3              : 2.14e-12 * exp(411/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H4     + O              -> HCO      + CH3              : 5.3e-12 * exp(-640/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H4     + H       + M    -> C2H5     + M                : k0 = 7.7e-30 * exp(-380/Tn) && \
                                                                                       kinf = 6.6e-15 * (Tn)^(1.28) * exp(-650/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H      + O2             -> CO       + HCO              : 4.0e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H      + H2             -> C2H2     + H                : 1.2e-11 * exp(-491/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H      + CH4            -> C2H2     + CH3              : 1.2e-11 * exp(-491/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H      + C2H6           -> C2H2     + C2H5             : 3.5e-11 * exp(2.9/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H      + H       + M    -> C2H2     + M                : k0 = 1.26e-18 * (Tn)^(-3.1) * exp(-721/Tn) && \
                                                                                       kinf = 3.0e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C        + OH             -> CO       + H                : 4.0e-11 * exp(2.9/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C        + H2      + M    -> ^3CH2    + M                : k0 = 7.0e-32 && \
                                                                                       kinf = 2.06e-11 * exp(-57/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C        + O2             -> CO       + O                : 3.3e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + H              -> C        + H2               : 1.3e-10 * exp(-80/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + O              -> CO       + H                : 6.6e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + H2             -> ^3CH2    + H                : 3.1e-10 * exp(-1650/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + H2      + M    -> CH3      + M                : k0 = 1.5e-23 * exp(-1650/Tn) && \
                                                                                       kinf = 8.55e-11 * (Tn)^(0.15) # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + O2             -> CO       + OH               : 2.75e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + CO2            -> HCO      + CO               : 5.9e-12 * exp(-350/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + CH4            -> C2H4     + H                : 3.96e-8 * (Tn)^(-1.04) * exp(-36.1/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + C2H2           -> C3H2     + H                : 1.59e-9 * (Tn)^(-0.233) * exp(-16/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + C2H4           -> CH3C2H   + H                : 3.87e-9 * (Tn)^(-0.546) * exp(-29.6/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + C2H4           -> CH2CCH2  + H                : 3.87e-9 * (Tn)^(-0.546) * exp(-29.6/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + O              -> CH       + OH               : 8.0e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + O              -> CO       + H     + H        : 1.20e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + H       + M    -> CH3      + M                : k0 = 3.4e-32 * exp(736/Tn) && \
                                                                                       kinf = 7.3e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + H              -> CH       + H2               : 3.54e-11 * (Tn)^(0.32) # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + CO      + M    -> CH2CO    + M                : k0 = 1.0e-28 && \
                                                                                       kinf = 1.5e-15 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + ^3CH2          -> C2H2     + H2               : 2.0e-11 * exp(-400/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + C2H2    + M    -> CH3C2H   + M                : k0 = 6.0e-29 * exp(1680/Tn) && \
                                                                                       kinf = 1.0e-11 * exp(-3330/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + C2H3           -> CH3      + C2H2             : 3.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + C2H5           -> CH3      + C2H4             : 3.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH2CO    + H              -> CH3      + CO               : 3.0e-11 * exp(-1700/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH2CO    + O              -> H2CO     + CO               : 1.3e-12 * exp(-680/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH2CCH2  + H       + M    -> CH3      + C2H2   + M       : k0 = 8.0e-24 * (Tn)^(-2) * exp(1225/Tn) && \
                                                                                       kinf = 9.7e-13 * exp(-1550/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH2CCH2  + H       + M    -> C3H5     + M                : k0 = 8.0e-24 * (Tn)^(-2) * exp(1225/Tn) && \
                                                                                       kinf = 6.6e-12 * exp(-1360/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + O2      + M    -> CH3O2    + M                : k0 = 4.0e-31 * (Tn/300)^(-3.6) && \
                                                                                       kinf = 1.2e-12 * (Tn/300)^(1.1) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + CO      + M    -> CH3CO    + M                : k0 = 1.26e-33 * exp(-1636/Tn) && \
                                                                                       kinf = 2.63e-13 * exp(-3007/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + OH             -> CO       + H2     + H2      : 6.7e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + HCO            -> CH4      + HCO              : 6.8e-12 * exp(-4450/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + C2H3           -> C3H5     + H                : 2.4e-13 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3O2    + H              -> CH4      + O2               : 1.4e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3O2    + H              -> H2O      + H2CO             : 1.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3O2    + O              -> H2CO     + HO2              : 1.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3CO    + H              -> CH4      + CO               : 1.0e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3CO    + O              -> H2CO     + HCO              : 5.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3CO    + CH3            -> C2H6     + CO               : 1.4e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3CO    + CH3            -> CH4      + CH2CO            : 1.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3CHO   + H              -> CH3CO    + H2               : 6.8e-15 * (Tn)^(1.16) * exp(-1210/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3CHO   + O              -> CH3CO    + OH               : 9.7e-12 * exp(-910/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3CHO   + OH             -> CH3CO    + H2O              : 1.95e-14 * (Tn)^(0.73) * exp(560/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3CHO   + CH3            -> CH3CO    + CH4              : 3.3e-30 * (Tn)^(5.64) * exp(-1240/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3C2H   + H       + M    -> CH3      + C2H2   + M       : k0 = 2.0e-29 && \
                                                                                       kinf = 3.98e-11 * exp(-1152/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3C2H   + H       + M    -> C3H5     + M                : k0 = 8.0e-24 * (Tn)^(2.0) * exp(-1225/Tn) && \
                                                                                       kinf = 6.0e-11 * exp(-1233/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2       + O              -> C        + CO               : 5.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C2       + O2             -> CO       + CO               : 1.5e-11 * exp(-550/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2       + H2             -> C2H      + H                : 1.77e-10 * exp(-1469/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2       + CH4            -> C2H      + CH3              : 5.05e-11 * exp(-297/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H      + O              -> CO       + CH               : 1.70e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H      + C3H8           -> C2H2     + C3H7             : 7.8e-11 * exp(3.0/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H2     + O              -> ^3CH2    + CO               : 3.0e-11 * exp(-1600/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H2     + OH      + M    -> C2H2OH   + M                : k0 = 5.5e-30 && \
                                                                                       kinf = 8.3e-13 * (Tn/300)^(2.0) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H2     + OH      + M    -> CH2CO    + H      + M       : k0 = 5.8e-31 * exp(-1258/Tn) && \
                                                                                       kinf = 1.4e-12 * exp(388/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H2OH   + H              -> H2O      + C2H2             : 5.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H2OH   + H              -> H2       + CH2CO            : 3.3e-11 * exp(-2000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H2OH   + O              -> OH       + CH2CO            : 3.3e-11 * exp(-2000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H2OH   + OH             -> H2O      + CH2CO            : 1.7e-11 * exp(-1000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H3     + O              -> CH2CO    + H                : 1.6e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H3     + OH             -> C2H2     + H2O              : 5.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H3     + CH3            -> C2H2     + CH4              : 3.4e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H3     + CH3     + M    -> C3H6     + M                : k0 = 6.0e-28 * exp(1680/Tn) && \
                                                                                       kinf = 1.2e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H3     + C2H3           -> C2H4     + C2H2             : 2.4e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H3     + C2H5           -> C2H4     + C2H4             : 8.0e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H3     + C2H5    + M    -> CH3      + C3H5   + M       : k0 = 1.9e-27 && \
                                                                                       kinf = 2.5e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H4     + OH      + M    -> C2H4OH   + M                : k0 = 1.0e-28 * (300/Tn)^(4.5) && \
                                                                                       kinf = 8.8e-12 * (300/Tn)^(0.85) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H4OH   + H              -> H2O      + C2H4             : 5.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H4OH   + H              -> H2       + CH3CHO           : 3.3e-11 * exp(-2000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H4OH   + O              -> OH       + CH3CHO           : 1.7e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H4OH   + OH             -> H2O      + CH3CHO           : 1.7e-11 * exp(-1000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + OH             -> CH3CHO   + H2               : 1.1e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + O              -> CH3CHO   + H                : 9.1e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + CH3            -> C2H4     + CH4              : 3.25e-11 * (Tn)^(-0.5) # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + C2H3           -> C2H6     + C2H2             : 8.0e-13 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + C2H5           -> C2H6     + C2H4             : 2.3e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + H       + M    -> C2H6     + M                : k0 = 5.5e-22 * (Tn)^(-2) * exp(-1040/Tn) && \
                                                                                       kinf = 1.66e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C2H5     + H              -> C2H4     + H2               : 3.0e-12 # Tian et al. [2011] ")
reaction_rate_list.append(" C3H2     + H       + M    -> C3H3     + M                : k0 = 2.52e-28 && \
                                                                                       kinf = 5.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C3H3     + H       + M    -> CH3C2H   + M                : k0 = 5.5e-27 && \
                                                                                       kinf = 1.15e-10 * exp(-276/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C3H3     + H       + M    -> CH2CCH2  + M                : k0 = 5.5e-27 && \
                                                                                       kinf = 2.5e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C3H5     + H              -> CH3C2H   + H2               : 1.4e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C3H5     + H       + M    -> C3H6     + M                : k0 = 2.0e-28 && \
                                                                                       kinf = 2.8e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C3H5     + H              -> CH4      + C2H2             : 1.5e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C3H5     + CH3            -> CH3C2H   + CH4              : 5.0e-12 * (Tn)^(-0.32) * exp(132/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C3H5     + CH3            -> CH2CCH2  + CH4              : 1.69e-10 * (Tn)^(-0.32) * exp(66/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C3H6     + OH             -> CH3CHO   + CH3              : 4.1e-12 * exp(540/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C3H6     + O              -> CH3      + CH3     + CO     : 4.1e-12 * exp(-38/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C3H6     + H       + M    -> C3H7     + M                : k0 = 1.3e-28 * exp(-380/Tn) && \
                                                                                       kinf = 2.2e-11 * exp(-785/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" C3H7     + CH3            -> C3H6     + CH4              : 1.9e-11 * (Tn)^(-0.32) # Tian et al. [2011] ")
reaction_rate_list.append(" C3H7     + OH             -> C2H5CHO  + H2               : 1.1e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" C3H7     + O              -> C2H5CHO  + H                : 1.1e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CH2CCH2  + H              -> CH3C2H   + H                : 1.1e-11 * exp(-1000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" O        + H2CO           -> OH       + HCO              : 3.4e-11 * exp(-1600/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + C2H2    + M    -> CH2CCH2  + M                : 1.0e-11 * exp(-3330/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" ^1CH2    + H2             -> ^3CH2    + H2               : 1.26e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" C3H5     + H              -> CH2CCH2  + H2               : 1.4e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" HCO      + H2CO           -> CH3O     + CO               : 1.0e-17 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3O     + CO             -> CH3      + CO2              : 2.6e-11 * exp(-5940/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HCO      + NO2            -> OH       + NO      + CO     : 5.1e-11*0.4 # Tian et al. [2011] ")
reaction_rate_list.append(" HCO      + NO2            -> NO       + H       + CO2    : 5.1e-11*0.6 # Tian et al. [2011] ")
reaction_rate_list.append(" CH4      + O(1D)          -> CH3O     + H                : 3.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" HCN      + OH             -> CN       + H2O              : 2.41e-11 * exp(-5503/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HCN      + OH             -> CO       + NH2              : 1.07e-13 * exp(-5893/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HCN      + OH      + M    -> HCNOH    + M                : k0 = 1.5e-31 * exp(-875/Tn) && \
                                                                                       kinf = 1.16e-13 * exp(-400/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HCN      + O              -> NH       + CO               : 8.88e-13 * (Tn/298)^(1.21) * exp(-3851/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HCN      + O              -> CN       + OH               : 1.43e-12 * (Tn/298)^(1.47) * exp(-3801/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + N              -> HCN      + H                : 2.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" ^3CH2    + O              -> CO       + H2               : 8.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH4      + N              -> HCN      + H2      + H      : 2.51e-14 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + N              -> HCN      + H       + H      : 3.32e-13 # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + N              -> HCN      + H2               : 4.8e-11 * exp(-420/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH3      + N              -> H2CN     + H                : 4.3e-10 * exp(-420/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H2CN     + NH2            -> HCN      + NH3              : 5.42e-11 * (300/Tn)^(1.06) * exp(-60.8/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H2CN     + O              -> HCN      + OH               : 8.3e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" H2CN     + H              -> HCN      + H2               : 8.3e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" H2CN     + N              -> HCN      + NH               : 1.0e-10 * exp(-200/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" H2CN     + M              -> HCN      + H       + M      : 5.5e-17 # Tian et al. [2011] ")
reaction_rate_list.append(" H2CN     + OH             -> HCN      + H2O              : 8.3e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" HCNOH    + NH2            -> HNCO     + NH3              : 3.3e-11 * exp(-2000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HCNOH    + O              -> HNCO     + OH               : 3.3e-11 * exp(-2000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HCNOH    + OH             -> HNCO     + H2O              : 1.7e-11 * exp(-1000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HCNOH    + H              -> HNCO     + H2               : 3.3e-11 * exp(-2000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HCNOH    + H              -> HCN      + H2O              : 5.0e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" HNCO     + OH             -> NCO      + H2O              : 0.9*6.1e-17 * (Tn)^(1.5) * exp(-1809/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HNCO     + OH             -> CO2      + NH2              : 0.1*6.1e-17 * (Tn)^(1.5) * exp(-1809/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HNCO     + H              -> NH2      + CO               : 8.63e-14 * (Tn/298)^(2.49) * exp(-1181/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HNCO     + H              -> NCO      + H2               : 2.67e-13 * (Tn/298)^(2.41) * exp(-6194/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HNCO     + CN             -> HCN      + NCO              : 2.51e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" HNCO     + NH2            -> NH3      + NCO              : 6.1e-12 * exp(-1203/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" HNCO     + O              -> CO2      + NH               : 5.01e-13 * (Tn/298)^(1.41) * exp(-4292/Tn) # Tian et al. [2011]")
reaction_rate_list.append(" HNCO     + O              -> NCO      + OH               : 6.08e-13 * (Tn/298)^(2.11) * exp(-5753/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" NCO      + O              -> CN       + O2               : 4.05e-10 * (Tn/298)^(-1.43) * exp(-3502/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" NCO      + O              -> CO       + NO               : 6.5e-11 * (Tn/298)^(-1.14) * exp(-3502/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" NCO      + N              -> CO       + N2               : 5.5e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" NCO      + H              -> CO       + NH               : 2.2e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" NCO      + H2             -> HNCO     + H                : 6.54e-14 * (Tn)^(2.58) * exp(-2700/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" NCO      + OH             -> HNCO     + O                : 5.38e-14 * (Tn/298)^(2.27) * exp(-497/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + NO             -> HCO      + O                : 1.37e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + NO             -> CO       + NH               : 2.0e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + NO             -> CN       + OH               : 1.4e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + NO             -> NCO      + H                : 4.4e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CH       + NO             -> HCO      + N                : 1.33e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + HNO            -> HCN      + NO               : 3.01e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + H       + M    -> HCN      + M                : k0 = 9.35e-30 * (Tn/298)^(-2.0) * exp(-521/Tn) && \
                                                                                       kinf = 1.73e-10 * (Tn/298)^(-0.5) # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + HNO2           -> HCN      + NO2              : 2.01e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + H2O            -> HCN      + OH               : 3.82e-11 * exp(-6704/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + NH3            -> HCN      + NH2              : 2.91e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + OH             -> NCO      + H                : 7.01e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + OH             -> HCN      + O                : 1.0e-11 * exp(-1000/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + HCO            -> NCO      + H                : 1.0e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + H2             -> HCN      + H                : 1.03e-11 * exp(-1771/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + C2H2           -> HCN      + C2H              : 2.19e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + C2H4           -> HCN      + C2H3             : 2.09e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + C2H6           -> HCN      + C2H5             : 2.08e-11 * (Tn/298)^(0.22) * exp(57.8/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + CH4            -> HCN      + CH3              : 5.11e-13 * (Tn/298)^(2.64) * exp(150/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + O              -> CO       + N                : 3.4e-11 * exp(-211/Tn) # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + N              -> N2       + C                : 3.01e-10 # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + O2             -> NCO      + O                : 2.3e-11 # Tian et al. [2011] ")
reaction_rate_list.append(" CN       + NO             -> CO       + N2               : 1.6e-13 # Tian et al. [2011] ")
reaction_rate_list.append(" NO       + C              -> CO       + N                : 3.49e-11 * (Tn/298)^(-0.02) # Tian et al. [2011] ")
reaction_rate_list.append(" NH       + N              -> N2       + H                : 1.95e-11 * (Tn/298)^(0.51) * exp(-9.63/Tn) # Tian et al. [2011] ")

##reaction_rate_list.append(" CO    + OH        -> CO2  + H        : 1.5e-13 # Chaffin et al. [2017] ")


#--------------------------------------------------------------------------------------------------------------------------
#
#                                                        Mars
#
#--------------------------------------------------------------------------------------------------------------------------
#Planet_list[planet_name, length of reaction_rate_list]
Planet_list.append(['Mars',len(reaction_rate_list)])
# special reaction ---------------------------------------------------------------------------------
reaction_rate_list.append("  CO2 + p* -> CO2+   + e- + p*      : Impact ionization ")
reaction_rate_list.append("  CO2 + e-* -> CO2+   + e- + e-*      : Impact ionization ")
reaction_rate_list.append("  N2  + p* -> N2+    + e- + p*      : Impact ionization ")
reaction_rate_list.append("  N2  + p* -> N+     + e- + p*      : Impact ionization ")
reaction_rate_list.append("  N2  + p* -> N           + p*      : Impact ionization ")
reaction_rate_list.append("  N2  + p* -> N(2D)       + p*      : Impact ionization ")
reaction_rate_list.append("  N2  + e* -> N2+    + e- + e*      : Impact ionization ")
reaction_rate_list.append("  N2  + e* -> N+     + e- + e*      : Impact ionization ")
reaction_rate_list.append("  N2  + e* -> N           + e*      : Impact ionization ")
reaction_rate_list.append("  N2  + e* -> N(2D)       + e*      : Impact ionization ")
reaction_rate_list.append("  CO2 -> CO2+   + e-   : 0.5e-17       # cosmic ray Molina-Cuberos et al. [2001] ")
reaction_rate_list.append("  CO2 -> O+(4S) + e-   : 0.1 * 0.5e-17 # cosmic ray Molina-Cuberos et al. [2001] ")
reaction_rate_list.append("  N2  -> N2+    + e-   : 0.5e-17       # cosmic ray Molina-Cuberos et al. [2001] ")
reaction_rate_list.append("  N2  -> N  + N(2D)    : 0.5e-17       # cosmic ray Molina-Cuberos et al. [2001] ")

# photoionization ---------------------------------------------------------------------------------
reaction_rate_list.append(" CO2 + hv -> CO2+   + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" CO2 + hv -> CO+    + e- + O  : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" CO2 + hv -> O+(4S)     + e- + CO : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" CO2 + hv -> C+     + e- + O2 : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" O2  + hv -> O2+    + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" O2  + hv -> O+(4S) + e- + O  : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" O   + hv -> O+(4S) + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" O   + hv -> O+(2D) + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" O   + hv -> O+(2P) + e-      : Photoionization @ Schunk and Nagy [2001] ")
#reaction_rate_list.append(" O   + hv -> O+(4P) + e-      : Photoionization @ Schunk and Nagy [2001] ") # NO LOSS for O+(4P)
reaction_rate_list.append(" CO  + hv -> CO+    + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" CO  + hv -> O+(4S)     + e- + C  : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" CO  + hv -> C+     + e- + O  : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" N2  + hv -> N2+    + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" N2  + hv -> N+     + e- + N(2D)  : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" H   + hv -> H+     + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" H2  + hv -> H2+    + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" H2  + hv -> H+     + e- + H  : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" He  + hv -> He+    + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" Ar  + hv -> Ar+    + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" H2O + hv -> H2O+   + e-      : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" H2O + hv -> OH+    + e- + H  : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" H2O + hv -> H+     + e- + OH : Photoionization @ Schunk and Nagy [2001] ")
reaction_rate_list.append(" H2O + hv -> O+(4S) + e- + H2 : Photoionization @ Schunk and Nagy [2001] ")

# photodissociation ---------------------------------------------------------------------------------
reaction_rate_list.append(" CO2    + hv -> CO     + O       : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" CO2    + hv -> CO     + O(1D)   : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" H2O    + hv -> H      + OH      : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" H2O    + hv -> H2     + O(1D)   : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" O3     + hv -> O2     + O       : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" O3     + hv -> O2     + O(1D)   : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" O2     + hv -> O      + O       : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" O2     + hv -> O      + O(1D)   : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" H2     + hv -> H      + H       : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" OH     + hv -> O      + H       : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" OH     + hv -> O(1D)  + H       : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" HO2    + hv -> OH     + O       : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" H2O2   + hv -> OH     + OH      : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" H2O2   + hv -> HO2    + H       : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" H2O2   + hv -> H2O    + O(1D)   : Photodissociation # Chaffin et al. [2017] ")
reaction_rate_list.append(" N2     + hv -> N      + N(2D)   : Photodissociation ")
reaction_rate_list.append(" NO     + hv -> N      + O       : Photodissociation ")
reaction_rate_list.append(" NO2    + hv -> NO     + O       : Photodissociation ")
reaction_rate_list.append(" NO2    + hv -> NO     + O(1D)   : Photodissociation ")
reaction_rate_list.append(" NO3    + hv -> NO     + O2      : Photodissociation ")
reaction_rate_list.append(" NO3    + hv -> NO2    + O       : Photodissociation ")
reaction_rate_list.append(" N2O    + hv -> N2     + O(1D)   : Photodissociation ")
reaction_rate_list.append(" N2O5   + hv -> NO3    + NO  + O : Photodissociation ")
reaction_rate_list.append(" N2O5   + hv -> NO3    + NO2     : Photodissociation ")
reaction_rate_list.append(" HONO   + hv -> OH     + NO      : Photodissociation ")
reaction_rate_list.append(" HNO3   + hv -> OH     + NO2     : Photodissociation ")
reaction_rate_list.append(" HNO3   + hv -> HONO   + O       : Photodissociation ")
reaction_rate_list.append(" HNO3   + hv -> HONO   + O(1D)   : Photodissociation ")
reaction_rate_list.append(" HO2NO2 + hv -> OH     + NO3     : Photodissociation ")
reaction_rate_list.append(" HO2NO2 + hv -> HO2    + NO2     : Photodissociation ")
reaction_rate_list.append(" H2CO   + hv -> HCO    + H       : Photodissociation ")
reaction_rate_list.append(" H2CO   + hv -> CO     + H2      : Photodissociation ")
reaction_rate_list.append(" H2CO   + hv -> CO     + H  + H  : Photodissociation ")

# Kim and Fox 1994 @ Jovian ionosphere : Recombination of Hn+, C+ with e- ---------------------------------------------------------------------------------
reaction_rate_list.append(" H+ + e- -> H : 4.0e-12* ( Te/250.0)^(-0.7) @ Bates and Dalgarno [1962] # RC1 in Kim and Fox [1994] ")
#reaction_rate_list.append(" He+ + e- -> He : 4.0e-12* ( Te/250.0)^(-0.7) @ Bates and Dalgarno [1962] # RC2 in Kim and Fox [1994] ")
#reaction_rate_list.append(" HeH+ + e- -> H +  He : 1.0e-8 * ( Te/300.0)^(-0.6) @ Yousif and Mitchell [1989] # RC3 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+ + e- -> H   +  H       : 2.3e-7 * ( Te/300.0)^(-0.4) @ Auerbach et al. [1977] # RC4 in Kim and Fox [1994] ")
reaction_rate_list.append(" H3+ + e- -> H2  +  H       : 4.4e-8 * ( Te/300.0)^(-0.5) @ Canosa et al. [1992] # RC5 in Kim and Fox [1994] ")
reaction_rate_list.append(" H3+ + e- -> H   +  H +  H  : 5.6e-8 * ( Te/300.0)^(-0.5) @ Mitchell et al. [1983] # RC6 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+ + H2 -> H3+ +  H       : 2.00e-9 @ Clow and Futrell [1972] # R1 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + e- -> C            : 4.0e-12* ( Te/250.0)^(-0.7) @ Bates and Dalgarno [1962] # RC7 in Kim and Fox [1994] ")

# Fox and Sung 2001 : ion - neutral ---------------------------------------------------------------------------------
reaction_rate_list.append(" CO2+   + O     -> CO     + O2+        : 1.64e-10 @ Fehsenfeld et al. [1970] # R1a in Fox and Sung [2001] ")
reaction_rate_list.append(" CO2+   + O     -> CO2    + O+(4S)         : 9.60e-11 @ Fehsenfeld et al. [1970] # R1b in Fox and Sung [2001] ")
reaction_rate_list.append(" CO2+   + O2    -> CO2    + O2+        : 5.50e-11 * (300/Ti)^0.82 for T = ~ 1500 [K] && \
                                                                    1.50e-11 * (1500/Ti)^(-0.75) for T = 1500 ~ [K] @ Anicich [1993a], Ferguson et al. [1992] # R2 in Fox and Sung [2001] ")
reaction_rate_list.append(" CO2+   + NO    -> NO+    + CO2        : 1.23e-10 @ Anicich [1993a] # R3 in Fox and Sung [2001] ")
reaction_rate_list.append(" CO2+   + N     -> NO     + CO+        : 3.40e-10 @ Scott et al. [1998] # R4 in Fox and Sung [2001] ")
reaction_rate_list.append(" CO2+   + N(2D) -> N+     + CO2        : 2.00e-10 @ estimated, see Fox [1982a] # R5 in Fox and Sung [2001] ")
reaction_rate_list.append(" CO2+   + H2    -> HCO2+  + H          : 8.70e-10 @ Scott et al. [1997] # R6 in Fox and Sung [2001] ")
reaction_rate_list.append(" CO2+   + H     -> HCO+   + O          : 4.46e-10 @ Scott et al. [1997] # R7a in Fox and Sung [2001] ")
reaction_rate_list.append(" CO2+   + H     -> H+     + CO2        : 2.35e-11 @ Scott et al. [1997] # R7b Fox and Sung [2001] ")
reaction_rate_list.append(" CO+    + O     -> CO     + O+(4S)         : 1.40e-10 @ Fehsenfeld and Ferguson [1972] # R8 in Fox and Sung [2001] ")
reaction_rate_list.append(" CO+    + NO    -> CO     + NO+        : 4.20e-10 @ Anicich [1993a] # R9 in Fox and Sung [2001] ")
reaction_rate_list.append(" CO+    + O2    -> O2+    + CO         : 1.50e-10 * (300/Ti)^1.1 @ Anicich [1993a], Miller et al. [1984b] # R10 in Fox and Sung [2001] ")
reaction_rate_list.append(" CO+    + CO2   -> CO2+   + CO         : 1.10e-9 @ Anicich [1993a] # R11 in Fox and Sung [2001] ")
reaction_rate_list.append(" CO+    + H2    -> HCO+   + H          : 7.50e-10 @ Scott et al. [1997] # R12a in Fox and Sung [2001] ")
#reaction_rate_list.append(" CO+    + H2    -> HOC+   + H          : 7.50e-10 @ Scott et al. [1997] # R12b in Fox and Sung [2001] ")  # NO LOSS REACTIONS FOR HOC+’
reaction_rate_list.append(" CO+    + H     -> H+     + CO         : 4.00e-10 @ Scott et al. [1997] # R13 in Fox and Sung [2001] ")
reaction_rate_list.append(" CO+    + N     -> NO+    + C          : 8.20e-11 @ Scott et al. [1998] # R14 in Fox and Sung [2001] ")
reaction_rate_list.append(" O2+    + N     -> NO+    + O          : 1.00e-10 @ Scott et al. [1998] # R15 in Fox and Sung [2001] ")
reaction_rate_list.append(" O2+    + N(2D) -> NO+    + O          : 1.80e-10 @ Goldan et al. [1966] # R16a in Fox and Sung [2001] ")
reaction_rate_list.append(" O2+    + N(2D) -> N+     + O2         : 8.65e-11 @ reverse of (R31b) O'Keefe et al. [1986] # R16b in Fox and Sung [2001] ")
reaction_rate_list.append(" O2+    + NO    -> NO+    + O2         : 4.50e-10 @ Midey and Viggiano [1999] # R17 in Fox and Sung [2001] ")

#reaction_rate_list.append(" O2+    + CO2    -> CO2+    + O2          : Is there a reaction like this？')
# very low CO2+ density at night: now it is avoided by setting NON ZERO ionizing source at night
# according to Girazian+ 2017, CO2+ at night should be larger than 1/cc

reaction_rate_list.append(" O2+    + C     -> CO+    + O          : 5.00e-11 @ Prasad and Huntress [1980], estimate  # R18a in Fox and Sung [2001] ")
reaction_rate_list.append(" O2+    + C     -> C+     + O2         : 5.00e-11 @  # R18b in Fox and Sung [2001] ")
reaction_rate_list.append(" O2+    + N2    -> NO+    + NO         : 1.00e-15 @ Ferguson [1973] # R19 in Fox and Sung [2001] ")
reaction_rate_list.append(" N2+    + N     -> N+     + N2         : 1.00e-11 @ Ferguson [1973] # R20 in Fox and Sung [2001] ")
reaction_rate_list.append(" N2+    + CO2   -> N2     + CO2+       : 9.00e-10 * (300/Ti)^0.23 @ Dotan et al. [2000] # R21 in Fox and Sung [2001] ")
reaction_rate_list.append(" N2+    + CO    -> N2     + CO+        : 7.60e-11 @ Frost et al. [1998] # R22 in Fox and Sung [2001] ")
reaction_rate_list.append(" N2+    + O2    -> N2     + O2+        : 5.10e-11 * (300/Ti)^1.16 for T = ~ 1000[K]&& \
                                                                    1.26e-11 * (1000/Ti)^(-0.67) for T = 1000 ~ 2000[K] && \
                                                                    2.39e-11 for T = 2000 ~ [K] @ Scott et al. [1999], Dotan et al. [1997] # R23 in Fox and Sung [2001] ")
reaction_rate_list.append(" N2+    + O     -> NO+    + N(2D)      : 1.33e-10 * (300/Ti)^0.44 for T = ~ 1500[K]  && \
                                                                    6.55e-11 * (1500/Ti)^(-0.2) for T = 1500 ~ [K] @ Scott et al. [1999], Mcfarland et al. [1974] # R24a in Fox and Sung [2001] ")
reaction_rate_list.append(" N2+    + O     -> O+(4S)     + N2         : 7.00e-12 * (300/Ti)^0.23 for T = ~ 1500[K] && \
                                                                    4.83e-12*(1500/Ti)^(-0.41) for T = 1500 ~ [K] @ Scott et al. [1999], Mcfarland et al. [1974] # R24b Fox and Sung [2001] ")
reaction_rate_list.append(" N2+    + NO    -> N2     + NO+        : 3.60e-10 @ Scott et al. [1999] # R25 in Fox and Sung [2001] ")
reaction_rate_list.append(" N2+    + Ar    -> Ar+    + N2         : 1.10e-11 * exp(-2089/Ti) @ reverse of (R71)  # R26 in Fox and Sung [2001] ")
reaction_rate_list.append(" N2+    + H2    -> N2H+   + H          : 1.52e-9 @ Uiterwaal et al. [1995] # R27 in Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + CO2   -> CO+    + NO         : 2.02e-10 @ Anicich [1993a]  # R28a in Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + CO2   -> N      + CO2+       : 9.18e-10 @ Anicich [1993a] # R28b in Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + NO    -> N      + NO+        : 4.72e-10 * (300/Ti)^0.24 @ Anicich [1993a], Fahey et al. [1981a] # R29a in Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + NO    -> N2+    + O          : 8.33e-11 * (300/Ti)^0.24 @ Anicich [1993a], Fahey et al. [1981a] # R29b in Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + CO    -> NO+    + C          : 6.16e-11 * (300/Ti)^0.5 @ Anicich [1993a], Miller et al. [1984b] # R30a in Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + CO    -> CO+    + N          : 4.93e-10 * (300/Ti)^0.5 @ Anicich [1993a], Miller et al. [1984b] # R30b in Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + CO    -> C+     + NO         : 5.60e-12 * (300/Ti)^0.5 @ Anicich [1993a], Miller et al. [1984b] # R30c in Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + O2    -> O2+    + N      : 2.02e-10 * (300/Ti)^(-0.45) for T = ~ 1000 [K] && \
                                                                    3.49e-10 for T = 1000 ~ [K] @ Dotan et al. [1997], O'Keefe et al. [1986] # R31a in Fox and Sung [2001] ") # N = N(4S)
reaction_rate_list.append(" N+     + O2    -> O2+    + N(2D)      : 8.65e-11 * (300/Ti)^(-0.45) for T = ~ 1000[K] && \
                                                                    1.49e-10 for T = 1000 ~ [K] @ Dotan et al. [1997], O'Keefe et al. [1986] # R31b in Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + O2    -> NO+    + O      : 4.32e-11 * (300/Ti)^(-0.45) for T = ~ 1000[K] && \
                                                                    7.47e-11 for T = 1000 ~ [K] @ Dotan et al. [1997], O'Keefe et al. [1986] # R31c in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" N+     + O2    -> NO+    + O(1D)      : 1.75e-10 * (300/Ti)^(-0.45) for T = ~ 1000[K] && \
                                                                    3.02e-10 for T = 1000 ~ [K] @ Dotan et al. [1997], O'Keefe et al. [1986] # R31d in Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + O2    -> O+(4S)     + NO         : 4.34e-11 * (300/Ti)^(-0.45) for T = ~ 1000[K] && \
                                                                        7.53e-11 for T = 1000 ~ [K] @ Dotan et al. [1997], O'Keefe et al. [1986] # R31e Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + O     -> N      + O+(4S)         : 2.20e-12 @ Constantinides et al. [1979], Bates [1989] # R32 in Fox and Sung [2001] ")
reaction_rate_list.append(" N+     + H2    -> NH+    + H          : 5.00e-10 @ Anicich [1993a] # R33 in Fox and Sung [2001] ")   # NH+ NO LOSS REACTIONS
reaction_rate_list.append(" O+(4S)     + O2    -> O      + O2+        : 1.60e-11 * (300/Ti)^0.52 for T = ~ 900[K] && \
                                                                        9.00e-12 * (900/Ti)^(-0.92) for T = 900 ~ [K] @ Hierl et al. [1997] # R34 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(4S)     + NO    -> NO+    + O          : 7.00e-13 * (300/Ti)^0.66 for T = ~ 300[K] && \
                                                                        7.00e-13 * (300/Ti)^(-0.87) for T = 300 ~ [K] @ Dotan and Viggiano [1999] # R35 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(4S)     + CO2   -> O2+    + CO         : 1.10e-9 for T = ~ 800[K] && \
                                                                        1.10e-9 * (800/Ti)^0.39 for T = 800 ~ [K] @ Anicich [1993a], Viggiano et al. [1992] # R36 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(4S)     + N2    -> NO+    + N          : 1.20e-12 * (300/Ti)^0.45 for T = ~ 1000 [K] && \
                                                                        7.00e-13 * (1000/Ti)^(-2.12) for T = 1000 ~ [K] @ Hierl et al. [1997] # R37 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(4S)     + N(2D) -> N+     + O          : 1.30e-10 @ Constantinides et al. [1979], Bates [1989] # R38 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(4S)     + C     -> C+     + O          : 1.00e-10 @ estimate; see text # R39 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(4S)     + H2    -> OH+    + H          : 1.35e-9 @ Li et al. [1997a] # R40 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(4S)     + H     -> H+     + O          : 6.40e-10 @ Anicich [1993a] # R41 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2D) + CO2   -> O2+    + CO         : 6.00e-11 @ Viggiano et al. [1990] # R42a in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2D) + CO2   -> CO2+   + O          : 1.00e-9 @ Viggiano et al. [1990] # R42b in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2D) + CO    -> CO+    + O          : 1.30e-9 @ Glosik et al. [1978] # R43 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2D) + O2    -> O2+    + O          : 7.00e-10 @ Johnsen and Biondi [1980] # R44 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2D) + NO    -> NO+    + O          : 1.20e-9 @ Glosik et al. [1978] # R45 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2D) + O     -> O+(4S) + O          : 1.00e-11 @ Torr and Torr [1980] # R46 in Fox and Sung [2001] ")##
reaction_rate_list.append(" O+(2D) + N2    -> N2+    + O          : 5.70e-10 * exp(-400/Ti) @ Li et al. [1997b] # R47 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2D) + N     -> N+     + O          : 1.50e-10 @ Dalgarno [1979] # R48 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2D) + H2    -> OH+    + H          : 1.50e-9 @ Li et al. [1997a] # R49a in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2D) + H2    -> H2+    + O          : 4.50e-11 @ Li et al. [1997a] # R49b in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2D) + H2    -> H+     + OH         : 1.50e-11 @ Li et al. [1997a] # R49c in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2D) + e-    -> O+(4S) + e-         : 6.03e-8 * (300/Te)^0.5 @ McLaughlin and Bell [1998] # R50 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + CO2   -> CO     + O2+        : 6.00e-11 @ Viggiano et al. [1990] # R51a in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + CO2   -> CO2+   + O          : 1.00e-9 @ Viggiano et al. [1990] # R51b in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + CO    -> CO+    + O          : 1.30e-9 @ Glosik et al. [1978] # R52 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + O2    -> O+(4S)     + O2         : 1.30e-10 @ Glosik et al. [1978] # R53a in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + O2    -> O2+    + O          : 1.30e-10 @ Glosik et al. [1978] # R53b in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + O     -> O+(2D) + O          : 5.20e-10 @ Rusch et al. [1977] # R54 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + N2    -> O+(4S)     + N2         : 6.20e-10 * exp(-340/Ti) @ Li et al. [1997b] # R55 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + N     -> O+(4S)     + N(2D)      : 1.00e-11 @ A. Dalgarno (private communication to Fox [1982a] # R56 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + NO    -> NO+    + O          : 1.20e-9 @ Glosik et al. [1978] # R57 in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + H2    -> OH+    + H          : 8.50e-10 @ Li et al. [1997a] # R58a in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + H2    -> H2+    + H          : 3.33e-10 @ Li et al. [1997a] # R58b in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + H2    -> H+     + OH         : 6.93e-11 @ Li et al. [1997a] # R58c in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + H2    -> H+     + O     + H  : 6.93e-11 @ Li et al. [1997a] # R58d in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + e-    -> O+(4S) + e-         : 3.03e-8 * (300/Te)^0.5 @ McLaughlin and Bell [1998] # R59a in Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) + e-    -> O+(2D) + e-         : 1.84e-7 * (300/Te)^0.5 @ McLaughlin and Bell [1998] # R59b in Fox and Sung [2001] ")
reaction_rate_list.append(" C+     + CO2   -> CO+    + CO         : 1.10e-9 @ Fahey et al. [1981b] # R60 in Fox and Sung [2001] ")
reaction_rate_list.append(" C+     + NO    -> NO+    + C          : 7.50e-10 * (300/Ti)^0.2 @ Anicich [1993a], Miller et al. [1984b] # R61 in Fox and Sung [2001] ")
reaction_rate_list.append(" C+     + O2    -> O+(4S)     + CO         : 5.22e-10 @ Anicich [1993a], Miller et al. [1984b] # R62a in Fox and Sung [2001] ")
reaction_rate_list.append(" C+     + O2    -> CO+    + O          : 3.48e-10 @ Anicich [1993a], Miller et al. [1984b] # R62b in Fox and Sung [2001] ")
reaction_rate_list.append(" C+     + H2    -> CH+    + H          : 7.40e-10 * exp(-4538/Ti) @ Hierl et al. [1997] # R63 in Fox and Sung [2001] ")
reaction_rate_list.append(" H+     + CO2   -> HCO+   + O          : 3.80e-9 @ Anicich [1993a] # R64 in Fox and Sung [2001] ")
reaction_rate_list.append(" H+     + O2    -> O2+    + H          : 1.17e-9 @ Huntress et al. [1974] # R65 in Fox and Sung [2001] ")
reaction_rate_list.append(" H+     + NO    -> H      + NO+        : 1.90e-9 @ Fehsenfeld and Ferguson [1972] # R66 in Fox and Sung [2001] ")
reaction_rate_list.append(" H+     + O     -> O+(4S)     + H          : 3.75e-10 ") # Special ???? after Schunk and nagy 2009
reaction_rate_list.append(" Ar+    + CO2   -> Ar     + CO2+       : 5.0e-10 for T = ~ 700[K] && \
                                                                    5.0e-10 * (700/Ti)^1 for T = 700 ~ [K] @ Dotan et al. [1999] # R68 in Fox and Sung [2001] ")
reaction_rate_list.append(" Ar+    + O2    -> Ar     + O2+        : 4.90e-11 * (300/Ti)^0.78 for T = ~ 900 [K] && \
                                                                    2.08e-11 * (900/Ti)^(-1.65) for T = 900 ~ [K] @ Midey and Viggiano [1998] # R69 in Fox and Sung [2001] ")
reaction_rate_list.append(" Ar+    + CO    -> Ar     + CO+        : 3.70e-11 * (300/Ti)^0.43 for T = ~ 900 [K] && \
                                                                    2.30e-11 * (900/Ti)^(-1) for T = 900 ~ [K] @ Midey and Viggiano [1998] # R70 in Fox and Sung [2001] ")
reaction_rate_list.append(" Ar+    + N2    -> Ar     + N2+        : 1.10e-11 * (300/Ti)^(-1.13) @ Anicich [1993a], Dotan and Lindinger [1982] # R71 in Fox and Sung [2001] ")
reaction_rate_list.append(" Ar+    + NO    -> Ar     + NO+        : 3.10e-10 @ Anicich [1993a] # R72 in Fox and Sung [2001] ")
reaction_rate_list.append(" Ar+    + H2    -> H2+    + Ar         : 1.78e-11 @ Anicich [1993a] # R73a in Fox and Sung [2001] ")
reaction_rate_list.append(" Ar+    + H2    -> ArH+   + H          : 8.72e-10 @ Anicich [1993a] # R73b in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + CO    -> C+     + O     + He : 1.60e-9 @ Anicich [1993a] # R74 in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + CO2   -> C+     + O2    + He : 2.00e-11 @ Anicich [1993a] # R75a in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + CO2   -> CO+    + O     + He : 7.80e-10 @ Anicich [1993a] # R75b in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + CO2   -> O+(4S)     + CO    + He : 1.40e-10 @ Anicich [1993a] # R75c in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + CO2   -> CO2+   + He         : 5.00e-11 @ Anicich [1993a] # R75d in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + O2    -> O+(2D) + O     + He : 2.37e-10 @ Gerlich [1992], Bischof and Linder [1986] # R76a in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + O2    -> O+(4S)     + O     + He : 2.39e-11 @ Gerlich [1992], Bischof and Linder [1986] # R76b in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + O2    -> O2+    + O          : 9.20e-12 @ Gerlich [1992], Bischof and Linder [1986] # R76c in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + O2    -> O+(2P) + O     + He : 6.04e-10 @ Gerlich [1992], Bischof and Linder [1986] # R76d in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + O2    -> O+(4S)     + O(1D) + He : 4.60e-11 @ Gerlich [1992], Bischof and Linder [1986] # R76e in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + O     -> O+(4S)     + He         : 1.00e-13 @ estimated; compare Dalgarno and Fox [1994] # R77 in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + N2    -> N+     + N     + He : 7.80e-10 @ Anicich [1993a] # R78a in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + N2    -> N2+    + He         : 5.20e-10 @ Anicich [1993a] # R78b in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + NO    -> N+     + O     + He : 1.35e-9 @ Anicich [1993a] # R79a in Fox and Sung [2001] ")
reaction_rate_list.append(" He+    + NO    -> O+(4S)     + N     + He : 1.00e-10 @ Anicich [1993a] # R79b in Fox and Sung [2001] ")

# Fox and Sung 2001 : neutral - neutral
#reaction_rate_list.append(" N     + CO2 -> NO    + CO      : 1.7e-16 @ see text # R80 in Fox and Sung [2001] ")
reaction_rate_list.append(" N     + O2  -> NO    + O       : 1.5e-14 * (1/Tn)^(-1) * exp(-3270/Tn) @ Baulch et al. [1994] # R81 in Fox and Sung [2001] ")
reaction_rate_list.append(" N     + NO  -> N2    + O       : 3.4e-11 @ Lee et al. [1978] # R82 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2D) + CO2 -> NO    + CO      : 3.60e-13 @ Herron [1999] # R83 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2D) + CO  -> N + CO      : 1.90e-12 @ Herron [1999] # R84 in Fox and Sung [2001] ") # N = N(4S)
reaction_rate_list.append(" N(2D) + O2  -> NO    + O(1D)   : 9.70e-12 * exp(-185/Tn) @ Herron [1999], Shihira et al. [1994]  # R85 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2D) + O   -> N + O       : 6.90e-13 @ Fell et al. [1990] # R86 in Fox and Sung [2001] ") # N = N(4S)
reaction_rate_list.append(" N(2D) + N2  -> N + N2      : 1.70e-14 @ Herron [1999] # R87 in Fox and Sung [2001] ") # N = N(4S)
reaction_rate_list.append(" N(2D) + NO  -> N2    + O       : 6.70e-11 @ Fell et al. [1990] # R88 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2D) + H2  -> NH    + H       : 4.20e-11 * exp(-880/Tn) @ Herron [1999] # R89 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2D) + e-  -> N + e-      : 3.86e-10 * (300/Te)^(-0.81) @ Berrington and Burke [1981] # R90 in Fox and Sung [2001] ") # N = N(4S)
reaction_rate_list.append(" N(2P) + CO2 -> N(2D) + CO2     : 2.00e-15 @ Herron [1999] # R91 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2P) + CO  -> N(2D) + CO      : 6.00e-15 @ Herron [1999] # R92 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2P) + O2  -> NO    + O   : 1.03e-12 * exp(-60/Tn) @ Shihira et al. [1994], branching ratios from Rawlins et al. [1989] # R93a in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" N(2P) + O2  -> NO    + O(1D)   : 1.03e-12 * exp(-60/Tn) @ Shihira et al. [1994], branching ratios from Rawlins et al. [1989] # R93b in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2P) + O2  -> NO    + O(1S)   : 1.03e-12 * exp(-60/Tn) @ Shihira et al. [1994], branching ratios from Rawlins et al. [1989] # R93c in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2P) + O   -> N(2D) + O       : 1.70e-11 @ Piper [1993] # R94 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2P) + NO  -> N(2D) + NO      : 2.90e-11 @ Herron [1999] # R95 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2P) + N2  -> N(2D) + N2      : 5.00e-17 @ Herron [1999] # R96 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2P) + N   -> N(2D) + N       : 6.20e-13 @ Young and Dunn [1975] # R97 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2P) + H2  -> N(2D) + H2      : 2.50e-15 @ Herron [1999] # R98 in Fox and Sung [2001] ")
reaction_rate_list.append(" N(2P) + e-  -> N + e-      : 2.04e-10 * (300/Te)^(-0.85) @ Berrington and Burke [1981] # R99a in Fox and Sung [2001] ") # N = N(4S)
reaction_rate_list.append(" N(2P) + e-  -> N(2D) + e-      : 9.50e-9 @ Berrington and Burke [1981] # R99b in Fox and Sung [2001] ")
reaction_rate_list.append(" C     + CO2 -> CO    + CO      : 7.62e-14 * (300/Tn)^(-0.5) * exp(-3480/Tn) @ McElroy and McConnell [1971], estimate # R100 in Fox and Sung [2001] ")
reaction_rate_list.append(" C     + NO  -> CN    + O       : 7.50e-11 * (300/Tn)^0.16 @ Chastaing et al. [2000] # R101a in Fox and Sung [2001] ")
reaction_rate_list.append(" C     + NO  -> CO    + N       : 7.50e-11 * (300/Tn)^0.16 @ Chastaing et al. [2000] # R101b in Fox and Sung [2001] ")
reaction_rate_list.append(" C     + O2  -> CO    + O       : 4.90e-11 * (300/Tn)^0.32 @ Chastaing et al. [2000] # R102 in Fox and Sung [2001] ")
reaction_rate_list.append(" O(1D) + CO2 -> O + CO2     : 6.80e-11 * exp(117/Tn) @ Streit et al. [1976] # R103 in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" O(1D) + CO  -> O + CO      : 3.60e-11 @ Schofield [1978] # R104 in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" O(1D) + O2  -> O + O2      : 3.20e-11 * exp(67/Tn) @ Atkinson et al. [1997] # R105 in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" O(1D) + O   -> O + O       : 6.47e-12 * (300/Tn)^(-0.14) @ Jamieson et al. [1992] # R106 in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" O(1D) + N2  -> O + N2      : 1.80e-11 * exp(107/Tn) @ Atkinson et al. [1997] # R107 in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" O(1D) + H2  -> OH    + H       : 1.10e-10 @ Atkinson et al. [1997] # R108 in Fox and Sung [2001] ")
reaction_rate_list.append(" O(1D) + e-  -> O + e-      : 2.87e-10 * (300/Te)^(-0.91) @ Berrington and Burke [1981] # R109 in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" O(1S) + CO2 -> O(1D) + CO2     : 2.02e-11 * exp(-1327/Tn) @ Capetanakis et al. [1993] # R110a in Fox and Sung [2001] ")
reaction_rate_list.append(" O(1S) + CO2 -> O + CO2     : 1.19e-11 * exp(-1327/Tn) @ Capetanakis et al. [1993] # R110b in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" O(1S) + O2  -> O(1D) + O2      : 1.36e-12 * exp(-815/Tn) @ Capetanakis et al. [1993] # R111a in Fox and Sung [2001] ")
reaction_rate_list.append(" O(1S) + O2  -> O + O2      : 3.04e-12 * exp(-815/Tn) @ Capetanakis et al. [1993] # R111b in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" O(1S) + O   -> O(1D) + O       : 0.00 @  # R112 in Fox and Sung [2001] ") # ?????????
reaction_rate_list.append(" O(1S) + N2  -> O(1D) + N2      : 5.00e-17 @ Atkinson and Welge [1972] # R113 in Fox and Sung [2001] ")
reaction_rate_list.append(" O(1S) + CO  -> O(1D) + CO      : 7.40e-14 * exp(-961/Tn) @ Capetanakis et al. [1993] # R114 in Fox and Sung [2001] ")
reaction_rate_list.append(" O(1S) + H2  -> O(1D) + H2      : 2.86e-16 @ Capetanakis et al. [1993] # R115 in Fox and Sung [2001] ")
reaction_rate_list.append(" O(1S) + e-  -> O(1D) + e-      : 8.50e-9 @ Berrington and Burke [1981] # R116a in Fox and Sung [2001] ")
reaction_rate_list.append(" O(1S) + e-  -> O + e-      : 1.56e-10 * (300/Te)^(-0.94) @ Berrington and Burke [1981] # R116b in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" H     + H   +  CO2 -> H2 + CO2 : 1.20e-32 * (300/Tn)^1.3 @ Tsang and Hampson [1986] # R117 in Fox and Sung [2001] ")

# Fox and Sung 2001 : Dissociative recombination
reaction_rate_list.append(" CO2+ + e- -> CO    + O     : 3.50e-7 * (300/Te)^0.5 @ Gougousi et al. [1997] # R118 in Fox and Sung [2001] ")
reaction_rate_list.append(" CO+  + e- -> C     + O     : 1.80e-7 * (300/Te)^0.55 @ Rosen et al. [1998] # R119a in Fox and Sung [2001] ")
reaction_rate_list.append(" CO+  + e- -> C     + O(1D) : 0.25e-7 * (300/Te)^0.55 @ Rosen et al. [1998] # R119b in Fox and Sung [2001] ")
reaction_rate_list.append(" CO+  + e- -> C(1D) + O     : 0.70e-7 * (300/Te)^0.55 @ Rosen et al. [1998] # R119c in Fox and Sung [2001] ")
reaction_rate_list.append(" O2+  + e- -> O + O : 0.39e-7 * (300/Te)^0.70 for T = ~ 1200 [K] && \
                                                         1.48e-8 * (1200/Te)^0.56 for T = 1200 ~ [K] @ Alge et al. [1983], Mehr and Biondi [1969], branching ratios from Kella et al. [1997] # R120a in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" O2+  + e- -> O + O(1D) : 0.86e-7 * (300/Te)^0.70 for T = ~ 1200 [K] && \
                                                         3.25e-8 * (1200/Te)^0.56 for T = 1200 ~ [K] @ Alge et al. [1983], Mehr and Biondi [1969], branching ratios from Kella et al. [1997] # R120b in Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" O2+  + e- -> O(1D) + O(1D) : 6.05e-8 * (300/Te)^0.70 for T = ~ 1200 [K] && \
                                                         2.29e-8 * (1200/Te)^0.56 for T = 1200 ~ [K] @ Alge et al. [1983], Mehr and Biondi [1969], branching ratios from Kella et al. [1997] # R120c in Fox and Sung [2001] ")
reaction_rate_list.append(" O2+  + e- -> O(1D) + O(1S) : 9.75e-9 * (300/Te)^0.70 for T = ~ 1200 [K] && \
                                                         3.69e-9 * (1200/Te)^0.56 for T = 1200 ~ [K] @ Alge et al. [1983], Mehr and Biondi [1969], branching ratios from Kella et al. [1997] # R120d in Fox and Sung [2001] ")
reaction_rate_list.append(" N2+  + e- -> N     + N(2D) : 1.01e-7 * (300/Te)^0.39 @ Zipf [1980], branching ratios from Kella et al. [1996] # R121a in Fox and Sung [2001] ")
reaction_rate_list.append(" N2+  + e- -> N(2D) + N(2D) : 1.01e-7 * (300/Te)^0.39 @ Zipf [1980], branching ratios from Kella et al. [1996] # R121b in Fox and Sung [2001] ")
reaction_rate_list.append(" N2+  + e- -> N     + N(2P) : 1.76e-8 * (300/Te)^0.39 @ Zipf [1980], branching ratios from Kella et al. [1996] # R121c Fox and Sung [2001] ")
reaction_rate_list.append(" NO+  + e- -> N(2D) + O     : 3.40e-7 * (300/Te)^0.5 @ Vejby-Christensen et al. [1998] # R122a in Fox and Sung [2001] ")
reaction_rate_list.append(" NO+  + e- -> N     + O     : 0.60e-7 * (300/Te)^0.5 @ Vejby-Christensen et al. [1998] # R122b in Fox and Sung [2001] ")

# Fox and Sung 2001 : transition 
reaction_rate_list.append(" O+(2D) -> O+(4S) + hv : 4.85e-5 @ Seaton and Osterbrock [1957] # Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) -> O+(2D) + hv : 1.71e-1 @ Seaton and Osterbrock [1957] # Fox and Sung [2001] ")
reaction_rate_list.append(" O+(2P) -> O+(4S) + hv : 4.8e-2 @ Seaton and Osterbrock [1957] # Fox and Sung [2001] ")
reaction_rate_list.append(" N(2D)  -> N  + hv : 1.07e-5 @ Wiese et al. [1966] # Fox and Sung [2001] ") # N = N(4S)
reaction_rate_list.append(" N(2P)  -> N(2D)  + hv : 7.9e-2 @ Wiese et al. [1966] # Fox and Sung [2001] ")
reaction_rate_list.append(" N(2P)  -> N  + hv : 5.0e-3 @ Wiese et al. [1966] # Fox and Sung [2001] ") # N = N(4S)
reaction_rate_list.append(" O(1D)  -> O  + hv : 9.3e-3 @ Froese-Fischer and Saha [1983] # Fox and Sung [2001] ") # O = O(3P)
reaction_rate_list.append(" O(1S)  -> O(1D)  + hv : 1.06 @ Kernahan and Pang [1975] # Fox and Sung [2001] ")
reaction_rate_list.append(" O(1S)  -> O  + hv : 4.5e-2 @ Kernahan and Pang [1975] # Fox and Sung [2001] ") # O = O(3P)

# Mukundan et al., 2020 for loss of minor ions HCO2+, HCO+, N2H+, HNO+, OH+ ---------------------------------------------------------------------------------
reaction_rate_list.append(" N2H+   + CO    -> HCO+     + N2     : 8.8e-10 @ Anicich and Huntress [1986] # Mukundan et al. [2020] ")
reaction_rate_list.append(" N2H+   + O     -> OH+      + N2     : 1.4e-10 @ Anicich and Huntress [1986] # Mukundan et al. [2020] ")
reaction_rate_list.append(" N2H+   + NO    -> HNO+     + N2     : 3.4e-10 @ Anicich [1993b] # Mukundan et al. [2020] ")
reaction_rate_list.append(" N2H+   + CO2   -> HCO2+    + N2     : 1.0e-9 @ Le-Teuff et al. [1999] # Mukundan et al. [2020] ")
reaction_rate_list.append(" OH+    + CO    -> HCO+     + O      : 3.55e-10 @ Anicich [1993b] # Mukundan et al. [2020] ")
reaction_rate_list.append(" OH+    + CO    -> CO+      + OH     : 3.55e-10 @ Anicich [1993b] # Mukundan et al. [2020] ")
reaction_rate_list.append(" OH+    + O     -> O2+      + H      : 7.1e-10 @ Le-Teuff et al. [1999] # Mukundan et al. [2020] ")
reaction_rate_list.append(" OH+    + N2    -> N2H+     + O      : 2.4e-10 @ Anicich [1993b] # Mukundan et al. [2020] ")
reaction_rate_list.append(" OH+    + NO    -> HNO+     + O      : 6.11e-10 @ Jones et al. [1981] # Mukundan et al. [2020] ")
reaction_rate_list.append(" OH+    + NO    -> NO+      + OH     : 3.59e-10 @ Fox [2015] # Mukundan et al. [2020] ")
reaction_rate_list.append(" OH+    + CO2   -> HCO2+    + O      : 1.4e-9 @ Le-Teuff et al. [1999] # Mukundan et al. [2020] ")
reaction_rate_list.append(" OH+    + O2    -> O2+      + OH     : 3.8e-10 @ Jones et al. [1981] # Mukundan et al. [2020] ")
reaction_rate_list.append(" HNO+   + CO    -> HCO+     + NO     : 8.6e-10 @ Fox [2015] # Mukundan et al. [2020] ")
reaction_rate_list.append(" HNO+   + NO    -> NO+      + HNO    : 7e-10 @ Fox [2015] # Mukundan et al. [2020] ")
reaction_rate_list.append(" HNO+   + CO2   -> HCO2+    + NO     : 9.45e-10 @ Fox [2015] # Mukundan et al. [2020] ")
reaction_rate_list.append(" HCO2+  + O     -> HCO+     + O2     : 1.0e-9 @ Le-Teuff et al. [1999] # Mukundan et al. [2020] ")
reaction_rate_list.append(" HCO2+  + CO    -> HCO+     + CO2    : 7.8e-10 @ Prasad and Huntress [1980] # Mukundan et al. [2020] ")
reaction_rate_list.append(" HCO2+  + N2    -> N2H+     + CO2    : 1.37e-9 @ Anicich and Huntress [1986] # Mukundan et al. [2020] ")
reaction_rate_list.append(" N2H+   + e-    -> N2       + H      : 2.325e-7 * (300/Te)^0.84 @ Fox [2015] # Mukundan et al. [2020] ")
reaction_rate_list.append(" OH+    + e-    -> O        + H      : 3.75e-8 * (300/Te)^0.5 @ Matta et al. [2013] # Mukundan et al. [2020] ")
reaction_rate_list.append(" HCO+   + e-    -> H        + CO     : 2.00e-7 * (300/Te)^1 @ Fox [2015] # Mukundan et al. [2020] ")
reaction_rate_list.append(" HNO+   + e-    -> H        + NO     : 3.00e-7 * (300/Te)^0.5 @ Fox [2015] # Mukundan et al. [2020] ")
reaction_rate_list.append(" HCO2+  + e-    -> H        + CO2    : 3.4e-7 * (300/Te)^0.5 @ Fox [2015] # Mukundan et al. [2020] ")

# Anicich 1993a ---------------------------------------------------------------------------------
reaction_rate_list.append(" H2+    + CO2   -> HCO2+    + H      : 2.35e-9 @ Anicich [1993a] ")
reaction_rate_list.append(" H2+    + CO    -> HCO+     + H      : 0.77 * 2.90e-9 @ Anicich [1993a] ")
reaction_rate_list.append(" H2+    + H2O   -> H2O+     + H2     : 0.53 * 7.30e-9 @ Anicich [1993a] ")
reaction_rate_list.append(" H2+    + H2O   -> H3O+     + H      : 0.47 * 7.30e-9 @ Anicich [1993a] ")
reaction_rate_list.append(" H2+    + O2    -> O2+      + H2     : 0.29 * 2.70e-9 @ Anicich [1993a] ")
reaction_rate_list.append(" H2+    + O2    -> HO2+     + H      : 0.71 * 2.70e-9 @ Anicich [1993a] ")
reaction_rate_list.append(" H3+    + O     -> OH+      + H2     : 0.5 * 8.00e-9 @ Anicich [1993a] ") # branch ratio is unknown
reaction_rate_list.append(" H3+    + O     -> H2O+     + H      : 0.5 * 8.00e-9 @ Anicich [1993a] ") # branch ratio is unknown
reaction_rate_list.append(" H3+    + H2O   -> H3O+     + H2     : 5.30e-9 @ Anicich [1993a] ")
reaction_rate_list.append(" H3+    + O2    -> HO2+     + H2     : 6.70e-10 @ Anicich [1993a] ")
reaction_rate_list.append(" H3+    + CO    -> HCO+     + H2     : 0.94 * 1.85e-9 @ Anicich [1993a] ")
reaction_rate_list.append(" H3+    + CO    -> HOC+     + H2     : 0.06 * 1.85e-9 @ Anicich [1993a] ")
reaction_rate_list.append(" H3+    + CO2   -> HCO2+    + H2     : 2.50e-9 @ Anicich [1993a] ")

reaction_rate_list.append(" H2O+   + H2    -> H3O+     + H      : 7.60e-10 @ Anicich [1993a] ")
reaction_rate_list.append(" O+(4S) + H2O   -> H2O+     + O      : 2.60e-9 @ Anicich [1993a] ")


## Shinagawa+ 1987 ---------------------------------------------------------------------------------
reaction_rate_list.append(" O+(4S) + CO2 -> O2+ + CO : 9.4e-10 @ Schunk and Nagy [1980] # Shinagawa et al. [1987]")
reaction_rate_list.append(" O+(4S) + N2 -> NO+ + N : 1.2e-12 * (300/Tn)^1 @ Schunk and Nagy [1980] # Shinagawa et al. [1987]")
reaction_rate_list.append(" O+(4S) + H -> H+ + O : 2.5e-11 * (1/Tn)^(-0.5) @ Schunk and Nagy [1980] # Shinagawa et al. [1987]")
reaction_rate_list.append(" CO2+ + O -> CO + O2+ : 1.64e-10 @ Schunk and Nagy [1980] # Shinagawa et al. [1987]")
reaction_rate_list.append(" CO2+ + O -> CO2 + O+(4S) : 9.6e-11 @ Schunk and Nagy [1980] # Shinagawa et al. [1987]")
reaction_rate_list.append(" H+ + O -> O+(4S) + H : 2.2e-11 * (1/Ti)^(-0.5) @ Schunk and Nagy [1980] # Shinagawa et al. [1987]")
reaction_rate_list.append(" H+ + CO2 -> HCO+ + O : 3.0e-9 @ Schunk and Nagy [1980] # Shinagawa et al. [1987]")
reaction_rate_list.append(" CO2+ + e- -> CO + O : 1.14e-4 * (1/Te)^1 @ Schunk and Nagy [1980] # Shinagawa et al. [1987]")
reaction_rate_list.append(" O2+ + e- -> O + O : 1.6e-7 * (300/Te)^0.55 @ Schunk and Nagy [1980] # Shinagawa et al. [1987]")

## Molina-Cuberos+ 2001 ---------------------------------------------------------------------------------
# positive ion reactions
#reaction_rate_list.append(" C+       +  CO2        -> CO+   + CO2      : 1.56e-9  @ Ikezoe et al. [1986]  # R04 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO+      +  CO2        -> CO2+  + CO       : 1.42e-9  @ Ikezoe et al. [1986]  # R05 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO2+     +  CO         -> CO+   + CO2      : 1.9e-12  @ Ikezoe et al. [1986]  # R06 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO2+     +  O2         -> O2+   + CO2      : 1.0e-10  @ Ikezoe et al. [1986]  # R08 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" N+       +  CO2        -> CO+   + NO       : 2.5e-10  @ Ikezoe et al. [1986]  # R20 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" N+       +  CO2        -> CO2+  + N        : 7.5e-10  @ Ikezoe et al. [1986]  # R21 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" N2+      +  CO         -> CO+   + N2       : 7.4e-11  @ Ikezoe et al. [1986]  # R22 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" N2+      +  CO2        -> CO2+  + N2       : 7.7e-10  @ Ikezoe et al. [1986]  # R23 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" N2+      +  O2         -> O2+   + N2       : 6.0e-11  @ Ikezoe et al. [1986]  # R24 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O+(4S)   +  CO2        -> O2+   + CO       : 1.0e-9  @ Ikezoe et al. [1986]  # R25 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O+(4S)   +  N2         -> NO+   + N        : 1.2e-12  @ Ikezoe et al. [1986]  # R26 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+      +  NO         -> NO    + O2       : 4.4e-10  @ Ikezoe et al. [1986]  # R29 in Molina-Cuberos et al. [2001] ")

# positive cluster ions
reaction_rate_list.append(" CO2+       +  CO2    + M      -> CO2+(CO2)  + M          : 2.5e-28  @ Ikezoe et al. [1986]  # R07 in Molina-Cuberos et al. [2001] ")
reaction_rate_list.append(" CO2+(CO2)  +  O2              -> O2+        + CO2  + CO2 : 1.53e-10  @ Ikezoe et al. [1986]  # R09 in Molina-Cuberos et al. [2001] ")
reaction_rate_list.append(" CO2+(CO2)  +  O2              -> O2+(CO2)   + CO2        : 2.7e-11  @ Ikezoe et al. [1986]  # R10 in Molina-Cuberos et al. [2001] ")
##reaction_rate_list.append(" H3O+       +  H2O    + M      -> H3O+H2O    + M          : 3.4e-27  @ Ikezoe et al. [1986]  # R11 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+OH     +  H2O             -> H3O+(H2O)  + OH         : 1.4e-9   @ Ikezoe et al. [1986]  # R12 in Molina-Cuberos et al. [2001] ")
##reaction_rate_list.append(" H3O+H2O    +  H2O    + M      -> H3O+(H2O)2 + M          : 2.3e-27  @ Ikezoe et al. [1986]  # R11 in Molina-Cuberos et al. [2001] ")
reaction_rate_list.append(" O2+        +  CO2    + M      -> O2+(CO2)   + M          : 1.7e-29  @ Ikezoe et al. [1986]  # R27 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+        +  H2O    + M      -> O2+H2O     + M          : 2.8e-28  @ Ikezoe et al. [1986]  # R28 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+        +  O2     + M      -> O2+O2      + M          : 1.0e-30  @ Turunen et al. [1996] # R30 in Molina-Cuberos et al. [2001] ")
reaction_rate_list.append(" O2+(CO2)   +  CO2             -> O2+        + CO2  + CO2 : 2.4e-13  @ Ikezoe et al. [1986]  # R31 in Molina-Cuberos et al. [2001] ")
reaction_rate_list.append(" O2+(CO2)   +  H2O             -> O2+(H2O)   + CO2        : 1.1e-9   @ Ikezoe et al. [1986]  # R33 in Molina-Cuberos et al. [2001] ")
##reaction_rate_list.append(" O2+CO2     +  CO2    + M      -> O2+(CO2)2  + M          : 1.1e-9   @ Ikezoe et al. [1986]  # R33 in Molina-Cuberos et al. [2001] ")
##reaction_rate_list.append(" O2+(CO2)2  +  H2O             -> products                : 2.3e-9   @ Ikezoe et al. [1986]  # R34 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+H2O     +  H2O             -> H3O+       + OH   + O2  : 2.04e-10 @ Ikezoe et al. [1986]  # R35 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+H2O     +  H2O             -> H3O+OH     + O2         : 9.96e-10 @ Ikezoe et al. [1986]  # R36 in Molina-Cuberos et al. [2001] ")
##reaction_rate_list.append(" O2+H2O     +  H2O    + M      -> O2+(H2O)2  + M          : 1.3e-27  @ Ikezoe et al. [1986]  # R37 in Molina-Cuberos et al. [2001] ")
##reaction_rate_list.append(" O2+(H2O)2  +  H2O             -> products                : 3.24e-10 @ Ikezoe et al. [1986]  # R38 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+O2      +  CO2             -> O2+CO2     + O2         : 4.0e-13  @ Ikezoe et al. [1986]  # R39 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+O2      +  H2O             -> O2+H2O     + O2         : 1.7e-9   @ Ikezoe et al. [1986]  # R40 in Molina-Cuberos et al. [2001] ")
#
#  # Photodissociation of positive ions
#reaction_rate_list.append(" O2+O2       -> O2+    + O2         : 0.13   @ Totmatsu [1990]  # R41 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+H2O      -> O2+    + H2O        : 0.18   @ Totmatsu [1990]  # R42 in Molina-Cuberos et al. [2001] ")
##reaction_rate_list.append(" O2+(H2O)2   -> O2+    + 2H2O       : 0.27   @ Totmatsu [1990]  # R43 in Molina-Cuberos et al. [2001] ")
#
## electron attachment to neutrals
reaction_rate_list.append(" e- +  O         -> O-   + hv        : 1.3e-15                                             @ Banks and Kockarts [1973]  # R44 in Molina-Cuberos et al. [2001] ")
reaction_rate_list.append(" e- +  O2  + M   -> O2-  + M  + hv   : 2.0e-31 * (Tn/300)^1 * exp(-600/Tn) @ Turunen et al. [1996]  # R45 in Molina-Cuberos et al. [2001] ")
reaction_rate_list.append(" e- +  O3        -> O-   + O2        : 9.1e-12 * (Tn/300)^(-1.46)                  @ Turunen et al. [1996]  # R46 in Molina-Cuberos et al. [2001] ")
#
## negative ion reactions
#reaction_rate_list.append(" O-         +  CO2    + M      -> CO3-        + M          : 1.1e-27  @ Ikezoe et al. [1986]  # R47 in Molina-Cuberos et al. [2001] ")
reaction_rate_list.append(" O2-        +  CO2    + M      -> CO4-        + M          : 1.3e-29  @ Ikezoe et al. [1986]  # R48 in Molina-Cuberos et al. [2001] ")
reaction_rate_list.append(" O3-        +  CO2             -> CO3-        + O2         : 5.5e-10  @ Ikezoe et al. [1986]  # R49 in Molina-Cuberos et al. [2001] ")
##reaction_rate_list.append(" CO3-       +  NO              -> NO2-        + CO2        : 1.1e-11  @ Ikezoe et al. [1986]  # R50 in Molina-Cuberos et al. [2001] ")
##reaction_rate_list.append(" CO3-       +  NO2             -> NO3-        + CO2        : 2.0e-10  @ Ikezoe et al. [1986]  # R51 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO3-       +  H2O    + M      -> CO3-H2O     + M          : 1.0e-28  @ Fritzenwallner and Kopp [1998]  # R52 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO3-H2O    +  M               -> CO3-        + H2O   + M  : 7.2e-4 * (Tn/300)^1 * exp(-7050/Tn)   @ Fritzenwallner and Kopp [1998]  # R53 in Molina-Cuberos et al. [2001] ")
#
#reaction_rate_list.append(" CO3-H2O    +  H2O    + M      -> CO3-(H2O)2  + M          : 1.0e-28  @ Fritzenwallner and Kopp [1998]  # R57 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO3-(H2O)2 +  M               -> CO3-H2O     + H2O   + M  : 6.5e-3 * (Tn/300)^1 * exp(-6800/Tn)   @ Fritzenwallner and Kopp [1998]  # R58 in Molina-Cuberos et al. [2001] ")
#
#reaction_rate_list.append(" CO4-       +  O               -> CO3-        + O2         : 1.4e-10  @ Ikezoe et al. [1986]  # R60 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO4-       +  O3              -> O3-         + CO2   + O2 : 1.3e-10  @ Ikezoe et al. [1986]  # R60 in Molina-Cuberos et al. [2001] ")
#
## electron Photodetachment and photodissociation of negative ions
#reaction_rate_list.append(" CO-         -> e-     + CO     : 0.01    @ Whitten et al. [1971] estimate  # R73 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO3-        -> e-     + CO3    : 0.0098  @ Turunen et al. [1996]  # R74 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO3-        -> O-     + CO2    : 0.076   @ Turunen et al. [1996]  # R75 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO3-H2O     -> CO3-   + H2O    : 0.27    @ Turunen et al. [1996]  # R76 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO3-(H2O)2  -> CO3-   + 2H2O   : 0.27    @ estimate  # R77 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" CO4-        -> O2-    + CO2    : 0.0028  @ Turunen et al. [1996]  # R78 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O-          -> e-     + O      : 0.62    @ Turunen et al. [1996]  # R82 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2-         -> e-     + O2     : 0.17    @ Turunen et al. [1996]  # R83 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O3-         -> e-     + O3     : 0.021   @ Turunen et al. [1996]  # R84 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O3-         -> O-     + O2     : 0.21    @ Turunen et al. [1996]  # R85 in Molina-Cuberos et al. [2001] ")
#
## ion - ion recombination
#reaction_rate_list.append(" O2+CO2     + CO3-(H2O)2  -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+H2O     + CO3-(H2O)2  -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)  + CO3-(H2O)2  -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)2 + CO3-(H2O)2  -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)3 + CO3-(H2O)2  -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)4 + CO3-(H2O)2  -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)5 + CO3-(H2O)2  -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)6 + CO3-(H2O)2  -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)7 + CO3-(H2O)2  -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+CO2     + CO3-H2O     -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+H2O     + CO3-H2O     -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)  + CO3-H2O     -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)2 + CO3-H2O     -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)3 + CO3-H2O     -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)4 + CO3-H2O     -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)5 + CO3-H2O     -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)6 + CO3-H2O     -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)7 + CO3-H2O     -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+CO2     + CO3-        -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" O2+H2O     + CO3-        -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)  + CO3-        -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)2 + CO3-        -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)3 + CO3-        -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)4 + CO3-        -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)5 + CO3-        -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)6 + CO3-        -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#reaction_rate_list.append(" H3O+(H2O)7 + CO3-        -> products : 6.0e-8 * (300/Tn)^0.5   @ Arijs [1992]  # R86 in Molina-Cuberos et al. [2001] ")
#
#
#
## Solomon et al. [1981], Pavlov et al. [2014] ---------------------------------------------------------------------------------
#reaction_rate_list.append(" H3O+OH     +  e-              -> H  + OH  + H2O          : 1.4e-6 * (300/Te)^0.66 @ R129 in Pavlov et al. [2014] # R9 in Solomon et al. [1981] ")
#
## Pavlov 2014
#reaction_rate_list.append(" H3O+       +  H2O    + M      -> H3O+(H2O)  + M          : 5.1e-27 * (300/Tn)^3.6 && \
#                                                                                       5.1e-27 * (300/Tn)^3.9 for T > 300 [K] @ Hamon et al. [2005]; Good et al. [1970a] # R94 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)  +  H2O    + M      -> H3O+(H2O)2 + M          : 2.3e-27 * (300/Tn)^7.5 @ Lau et al. [1982]; Good et al. [1970a] # R96 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)2 +  H2O    + M      -> H3O+(H2O)3 + M          : 2.8e-27 * (300/Tn)^3.1 && \
#                                                                                       2.4e-27 * (300/Tn)^8.1 for T > 290 [K] @ Lau et al. [1982] # R98 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)3 +  H2O    + M      -> H3O+(H2O)4 + M          : 4.6e-27 * (300/Tn)^5.1 && \
#                                                                                       9.0e-28 * (300/Tn)^14 for T > 250 [K] @ Lau et al. [1982]; Good et al. [1970b] # R100 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)4 +  H2O    + M      -> H3O+(H2O)5 + M          : 1.0e-27 * (300/Tn)^7.8 && \
#                                                                                       1.2e-28 * (300/Tn)^15.3 for T > 225 [K] @ Lau et al. [1982]; Good et al. [1970b] # R102 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)5 +  H2O    + M      -> H3O+(H2O)6 + M          : 1.0e-27 * (300/Tn)^7.8 && \
#                                                                                       1.2e-28 * (300/Tn)^15.3 for T > 225 [K]  # R104 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)6 +  H2O    + M      -> H3O+(H2O)7 + M          : 1.0e-27 * (300/Tn)^7.8 && \
#                                                                                       1.2e-28 * (300/Tn)^15.3 for T > 225 [K]  # R106 in Pavlov [2014] ")
#
#reaction_rate_list.append(" H3O+(H2O)  + M                -> H3O+       + H2O  + M   : 0.074 * (300/Tn)^4.6 * exp(-16164/Tn) && \
#                                                                                       0.074 * (300/Tn)^4.9 * exp(-16164/Tn) for T > 300 [K]  # R95 in Pavlov [2014] (reverse reaction) ")
#reaction_rate_list.append(" H3O+(H2O)2 + M                -> H3O+(H2O)  + H2O  + M   : 2.6e-3 * (300/Tn)^8.5 * exp(-10272/Tn) # R97 in Pavlov [2014] (reverse reaction) ")
#reaction_rate_list.append(" H3O+(H2O)3 + M                -> H3O+(H2O)2 + H2O  + M   : 0.083 * (300/Tn)^4.1 * exp(-8661/Tn) && \
#                                                                                       0.071 * (300/Tn)^9.1 * exp(-8661/Tn) for T > 290 [K]  # R99 in Pavlov [2014] (reverse reaction) ")
#reaction_rate_list.append(" H3O+(H2O)4 + M                -> H3O+(H2O)3 + H2O  + M   : 0.020 * (300/Tn)^6.1 * exp(-6143/Tn) && \
#                                                                                       3.9e-3 * (300/Tn)^15 * exp(-6143/Tn) for T > 250 [K]  # R101 in Pavlov [2014] (reverse reaction) ")
#reaction_rate_list.append(" H3O+(H2O)5 + M                -> H3O+(H2O)4 + H2O  + M   : 5.7e-3 * (300/Tn)^8.8 * exp(-5841/Tn) && \
#                                                                                       6.7e-4 * (300/Tn)^16.3 * exp(-5841/Tn) for T > 225 [K]  # R103 in Pavlov [2014] (reverse reaction) ")
#reaction_rate_list.append(" H3O+(H2O)6 + M                -> H3O+(H2O)5 + H2O  + M   : 0.029 * (300/Tn)^8.8 * exp(-5640/Tn) && \
#                                                                                       3.4e-3 * (300/Tn)^16.3 * exp(-5640/Tn) for T > 225 [K]  # R105 in Pavlov [2014] (reverse reaction) ")
#reaction_rate_list.append(" H3O+(H2O)7 + M                -> H3O+(H2O)6 + H2O  + M   : 0.030 * (300/Tn)^8.8 * exp(-5187/Tn) && \
#                                                                                       3.5e-3 * (300/Tn)^16.3 * exp(-5187/Tn) for T > 225 [K]  # R107 in Pavlov [2014] (reverse reaction) ")
#
#reaction_rate_list.append(" H3O+       +  e-              -> H          + H2O        : 7.6e-7 * (300/Te)^0.83 && \
#                                                                                       1.05e-6 * (300/Te)^1.1 for T > 1000 [K] @ Neau et al. [2000] # R111 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)  +  e-              -> H          + 2H2O       : 1.4e-6 * (300/Te)^0.66 && \
#                                                                                       1.54e-6 * (300/Te)^0.85 for T > 500 [K] @ Ojekull et al. [2007] # R112 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)2 +  e-              -> H          + 3H2O       : 2.48e-6 * (300/Te)^0.763 @ Ojekull et al. [2007] # R113 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)3 +  e-              -> H          + 4H2O       : 5.5e-7 * (300/Te)^0.781 @ Ojekull et al. [2008] # R114 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)4 +  e-              -> H          + 5H2O       : 3.80e-6 * (300/Te)^0.687 @ Ojekull et al. [2008] # R115 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)5 +  e-              -> H          + 6H2O       : 2.25e-6 * (300/Te)^0.652 @ Ojekull et al. [2008] # R116 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)6 +  e-              -> H          + 7H2O       : 1.3e-6 * (300/Te)^0.652  # R116 in Pavlov [2014] ")
#reaction_rate_list.append(" H3O+(H2O)7 +  e-              -> H          + 8H2O       : 7.8e-7 * (300/Te)^0.652  # R116 in Pavlov [2014] ")

# N4+, NO+(H2O) reactions
reaction_rate_list.append(" N2+        +  N2     + M      -> N4+       + M           : 6.8e-29 * (300/Tn)^2.23 * (1-0.00824*(300/Tn)^0.89) @ Troe [2005] # R31 in Pavlov [2014]")
reaction_rate_list.append(" N4+        +  M               -> N2+       + N2     + M  : 6.5e-5*(300/Tn)^3.23*(1-0.00824*(300/Tn)^0.89)*exp(-1300/Tn) @ Reverse reaction # R32 in Pavlov [2014]")
reaction_rate_list.append(" N4+        +  O2              -> O2+       + 2N2         : 1.5e-10 @ Schulz and Armentrout [1991a] # R33 in Pavlov [2014]")
reaction_rate_list.append(" N4+        +  CO2             -> CO2+      + 2N2         : 7.0e-10 @ Smith et al. [1978b] # R34 in Pavlov [2014]")
reaction_rate_list.append(" N4+        +  H2O             -> H2O+      + 2N2         : 3.0e-10 @ Smith et al. [1978b] # R35 in Pavlov [2014]")
reaction_rate_list.append(" NO+(H2O)   +  M               -> NO+       + H2O    + M  : 3.5e-4 * (300/Tn)^3.837 * exp(-9316/Tn) @ Reverse reaction  # R54 in Pavlov [2014] ") # N2 -> M 
reaction_rate_list.append(" N4+        +  e-              -> 2N2                     : 2.6e-6 * (300/Tn)^5 * exp(-3872/Tn) @ Cao and Johnsen [1991] # R109 in Verronen et al. [2016] RNW")


# Verronen+ 2016: photochemistry in D region A2.1: positive ion reaction ---------------------------------------------------------------------------------
reaction_rate_list.append(" O2+        +  O2     + M      -> O4+         + M               : 4.0e-30 * (300/Tn)^2.93 @ Bohringer and Arnold [1982] # R6 in Verronen et al. [2016]")
reaction_rate_list.append(" O2+        +  H2O    + M      -> O2+(H2O)    + M               : 2.8e-28   # R7  in Verronen et al. [2016]")
reaction_rate_list.append(" O4+        +  H2O             -> O2+(H2O)    + O2              : 1.7e-9    # R10 in Verronen et al. [2016]")
reaction_rate_list.append(" O4+        +  O               -> O2+         + O3              : 3.0e-10   # R11 in Verronen et al. [2016]")
reaction_rate_list.append(" O2+(H2O)   +  H2O             -> H3O+(OH)    + O2              : 9.0e-10   # R32 in Verronen et al. [2016]")
reaction_rate_list.append(" O2+(H2O)   +  H2O             -> H+(H2O)     + OH      + O2    : 2.4e-10   # R33 in Verronen et al. [2016]")
reaction_rate_list.append(" H3O+(OH)   +  H2O             -> H+(H2O)2    + OH              : 2.0e-9    # R34 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)    +  H2O    + M      -> H+(H2O)2    + M               : 4.6e-27 * (300/Tn)^4    # R35 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)2   +  M               -> H+(H2O)     + H2O     + M     : 2.5e-2 * (300/Tn)^5 * exp(-15900/Tn)  # R36 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)2   +  H2O    + M      -> H+(H2O)3    + M               : 2.3e-27 * (300/Tn)^7.5    # R37 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)3   +  M               -> H+(H2O)2    + H2O     + M     : 2.6e-3 * (300/Tn)^8.5 * exp(-10272/Tn)  # R38 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)3   +  H2O    + M      -> H+(H2O)4    + M               : 3.6e-27 * (300/Tn)^8.1    # R39 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)4   +  M               -> H+(H2O)3    + H2O     + M     : 1.5e-1 * (300/Tn)^9.1 * exp(-9000/Tn)  # R40 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)4   +  H2O    + M      -> H+(H2O)5    + M               : 4.6e-28 * (300/Tn)^14    # R41 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)5   +  M               -> H+(H2O)4    + H2O     + M     : 1.7e-3 * (300/Tn)^15 * exp(-6400/Tn)  # R42 in Verronen et al. [2016]")

# Verronen+ 2016: photochemistry in D region A2.2: recombination of positive ions with e-
reaction_rate_list.append(" O4+        +  e-              -> 2O2           : 4.2e-6 * (300/Tn)^0.5   # R2 in Verronen et al. [2016]")
reaction_rate_list.append(" O2+(H2O)   +  e-              -> O2     + H2O  : 2.0e-6     # R10 in Verronen et al. [2016]")
reaction_rate_list.append(" H3O+(OH)   +  e-              -> OH + H + H2O  : 1.5e-6     # R11 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)    +  e-              -> H      + H2O  : 6.3e-7 * (300/Tn)^0.5    # R12 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)2   +  e-              -> H      + 2H2O : 2.5e-6 * (300/Tn)^0.1    # R13 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)3   +  e-              -> H      + 3H2O : 2.48e-6 * (300/Tn)^0.76    # R14 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)4   +  e-              -> H      + 4H2O : 3.6e-6     # R15 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)5   +  e-              -> H      + 5H2O : 5.0e-6     # R15 in Verronen et al. [2016]")

# Verronen+ 2016: photochemistry in D region A2.3: positive ion photodissociation
reaction_rate_list.append(" O2+(H2O)   -> O2+  + H2O      : 0.42   # R1 in Verronen et al. [2016]")

# Verronen+ 2016: photochemistry in D region A2.5: negative ion reactions
reaction_rate_list.append(" O-         +  O3              -> O3-         + O               : 8.0e-10   # R1 in Verronen et al. [2016]")
reaction_rate_list.append(" O-         +  CO2    + M      -> CO3-        + M               : 2.0e-28   # R5 in Verronen et al. [2016]")
reaction_rate_list.append(" O2-        +  O               -> O-          + O2              : 1.5e-10   # R11 in Verronen et al. [2016]")
reaction_rate_list.append(" O2-        +  O3              -> O3-         + O2              : 7.8e-10   # R12 in Verronen et al. [2016]")
#reaction_rate_list.append(" O2-        +  CO2    + M      -> CO4-        + M               : 3.2e-11   # R6 in Verronen et al. [2016]")
reaction_rate_list.append(" O2-        +  O2     + M      -> O4-         + M               : 3.4e-31   # R15 in Verronen et al. [2016]")
reaction_rate_list.append(" O3-        +  O               -> O2-         + O2              : 2.5e-10   # R20 in Verronen et al. [2016]")
reaction_rate_list.append(" O3-        +  CO2             -> CO3-        + O2              : 5.5e-10   # R22 in Verronen et al. [2016]")
reaction_rate_list.append(" O4-        +  O               -> O3-         + O2              : 4.0e-10   # R27 in Verronen et al. [2016]")
reaction_rate_list.append(" O4-        +  CO2             -> CO4-        + O2              : 4.3e-10   # R28 in Verronen et al. [2016]")
reaction_rate_list.append(" CO3-       +  O               -> O2-         + CO2             : 1.1e-10   # R35 in Verronen et al. [2016]")
reaction_rate_list.append(" CO3-       +  O2              -> O3-         + CO2             : 6.0e-15   # R36 in Verronen et al. [2016]")
reaction_rate_list.append(" CO3-       +  H2O    + M      -> CO3-(H2O)   + M               : 1.0e-28   # R40 in Verronen et al. [2016]")
reaction_rate_list.append(" CO4-       +  O3              -> O3-         + O2      + CO2   : 1.3e-10   # R40 in Verronen et al. [2016]")
reaction_rate_list.append(" CO4-       +  O               -> CO3-        + O2      + CO2   : 1.4e-10   # R40 in Verronen et al. [2016]")
reaction_rate_list.append(" CO3-(H2O)  +  M               -> CO3-        + H2O     + M     : 7.2e-4 * (300/Tn) * exp(-7050/Tn)   # R77 in Verronen et al. [2016]")
reaction_rate_list.append(" CO3-(H2O)  +  H2O    + M      -> CO3-(H2O)2  + M               : 1.0e-28   # R79 in Verronen et al. [2016]")
reaction_rate_list.append(" CO3-(H2O)2 +  M               -> CO3-(H2O)   + H2O     + M     : 6.5e-3 * (300/Tn) * exp(-6800/Tn)   # R81 in Verronen et al. [2016]")

# Verronen+ 2016: photochemistry in D region A2.6: electron detachment from negative ions
reaction_rate_list.append(" O-     +  O       -> O2     + e-      : 1.9e-10   # R1 in Verronen et al. [2016]")
reaction_rate_list.append(" O2-    +  O       -> O3     + e-      : 1.5e-10   # R6 in Verronen et al. [2016]")
reaction_rate_list.append(" O3-    +  O       -> 2O2    + e-      : 1.0e-10   # R10 in Verronen et al. [2016]")
reaction_rate_list.append(" O3-    +  O3      -> 3O2    + e-      : 1.0e-10   # R11 in Verronen et al. [2016]")

# Verronen+ 2016: photochemistry in D region A2.7: Photodetachment of electrons from negative ions
reaction_rate_list.append(" O-       -> O    + e-      : 1.4    # R1 in Verronen et al. [2016]")
reaction_rate_list.append(" O2-      -> O2   + e-      : 0.38   # R2 in Verronen et al. [2016]")
reaction_rate_list.append(" O3-      -> O3   + e-      : 4.7e-2 # R3 in Verronen et al. [2016]")

# Verronen+ 2016: photochemistry in D region A2.8: negative ion photodissociation
reaction_rate_list.append(" O3-       -> O-    + O2     : 0.47   # R1 in Verronen et al. [2016]")
reaction_rate_list.append(" O4-       -> O2-   + O2     : 0.24   # R2 in Verronen et al. [2016]")
reaction_rate_list.append(" CO3-      -> O-    + CO2    : 0.15   # R3 in Verronen et al. [2016]")
reaction_rate_list.append(" CO4-      -> O2-   + CO2    : 6.2e-3 # R4 in Verronen et al. [2016]")
reaction_rate_list.append(" CO3-(H2O) -> CO3-  + H2O    : 0.6    # R5 in Verronen et al. [2016]")

# Verronen+ 2016: photochemistry in D region A2.9: ion-ion recombination two body reaction
reaction_rate_list.append(" H+(H2O)4   +  CO3-       -> H     + 4H2O   + O   + CO2    : 6.0e-8 * (300/Tn)^0.5   # R2 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)4   +  O2-        -> H     + 4H2O   + O2           : 6.0e-8 * (300/Tn)^0.5   # R6 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)4   +  CO4-       -> H     + 4H2O   + O2  + CO2    : 6.0e-8 * (300/Tn)^0.5   # R7 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)4   +  CO3-(H2O)2 -> H     + 6H2O   + O   + CO2    : 6.0e-8 * (300/Tn)^0.5   # R9 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)4   +  CO3-(H2O)  -> H     + 5H2O   + O   + CO2    : 6.0e-8 * (300/Tn)^0.5   # R11 in Verronen et al. [2016]")

reaction_rate_list.append(" H+(H2O)5   +  CO3-       -> H     + 5H2O   + O   + CO2    : 6.0e-8 * (300/Tn)^0.5   # R18 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)5   +  O2-        -> H     + 5H2O   + O2           : 6.0e-8 * (300/Tn)^0.5   # R22 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)5   +  CO4-       -> H     + 5H2O   + O2  + CO2    : 6.0e-8 * (300/Tn)^0.5   # R23 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)5   +  CO3-(H2O)2 -> H     + 7H2O   + O   + CO2    : 6.0e-8 * (300/Tn)^0.5   # R25 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)5   +  CO3-(H2O)  -> H     + 6H2O   + O   + CO2    : 6.0e-8 * (300/Tn)^0.5   # R27 in Verronen et al. [2016]")

reaction_rate_list.append(" H+(H2O)3   +  CO3-       -> H     + 3H2O   + O   + CO2    : 6.0e-8 * (300/Tn)^0.5   # R32 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)3   +  O2-        -> H     + 3H2O   + O2           : 6.0e-8 * (300/Tn)^0.5   # R38 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)3   +  CO4-       -> H     + 3H2O   + O2  + CO2    : 6.0e-8 * (300/Tn)^0.5   # R39 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)3   +  CO3-(H2O)2 -> H     + 5H2O   + O   + CO2    : 6.0e-8 * (300/Tn)^0.5   # R41 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)3   +  CO3-(H2O)  -> H     + 4H2O   + O   + CO2    : 6.0e-8 * (300/Tn)^0.5   # R43 in Verronen et al. [2016]")

reaction_rate_list.append(" O2+        +  CO3-       -> O2    + O   + CO2             : 6.0e-8 * (300/Tn)^0.5   # R98 in Verronen et al. [2016]")
reaction_rate_list.append(" O2+        +  O2-        -> 2O2                           : 6.0e-8 * (300/Tn)^0.5   # R102 in Verronen et al. [2016]")
reaction_rate_list.append(" O2+        +  CO4-       -> 2O2   + CO2                   : 6.0e-8 * (300/Tn)^0.5   # R103 in Verronen et al. [2016]")
reaction_rate_list.append(" O2+        +  CO3-(H2O)2 -> O2    + O   + 2H2O   + CO2    : 6.0e-8 * (300/Tn)^0.5   # R105 in Verronen et al. [2016]")
reaction_rate_list.append(" O2+        +  CO3-(H2O)  -> O2    + O   + H2O    + CO2    : 6.0e-8 * (300/Tn)^0.5   # R107 in Verronen et al. [2016]")

# Verronen+ 2016: photochemistry in D region A2.10: ion-ion recombination three body reaction
reaction_rate_list.append(" H+(H2O)4   +  CO3-       + M   -> H     + 4H2O   + O   + CO2  + M   : 1.25e-25 * (300/Tn)^4   # R1 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)5   +  CO3-       + M   -> H     + 5H2O   + O   + CO2  + M   : 1.25e-25 * (300/Tn)^4   # R3 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)4   +  CO3-(H2O)2 + M   -> H     + 6H2O   + O   + CO2  + M   : 1.25e-25 * (300/Tn)^4   # R9 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)5   +  CO3-(H2O)2 + M   -> H     + 7H2O   + O   + CO2  + M   : 1.25e-25 * (300/Tn)^4   # R10 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)4   +  CO3-(H2O)  + M   -> H     + 5H2O   + O   + CO2  + M   : 1.25e-25 * (300/Tn)^4   # R11 in Verronen et al. [2016]")
reaction_rate_list.append(" H+(H2O)5   +  CO3-(H2O)  + M   -> H     + 6H2O   + O   + CO2  + M   : 1.25e-25 * (300/Tn)^4   # R12 in Verronen et al. [2016]")


# Verronen+ 2016: photochemistry in D region A2.1: positive ion reaction including RNW---------------------------------------------------------------------------------
reaction_rate_list.append(" NO+        +  N2         + M   -> NO+(N2)     + M              : 3.0e-31 * (300/Tn)^4.3                # R18 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  CO2        + M   -> NO+(CO2)    + M              : 1.4e-29 * (300/Tn)^4                  # R19 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  H2O        + M   -> NO+(H2O)    + M              : 1.35e-28 * (300/Tn)^2.83              # R20 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(N2)    +  CO2              -> NO+(CO2)    + N2             : 1.0e-9                                # R21 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(N2)    +  H2O              -> NO+(H2O)    + M              : 1.0e-9                                # R22 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(N2)    +  M                -> NO+         + N2     + M     : 1.5e-8 * (300/Tn)^4.3 * exp(-2093/Tn) # R23 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(CO2)   +  H2O              -> NO+(H2O)    + CO2            : 1.0e-9                                # R24 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(CO2)   +  M                -> NO+         + CO2    + M     : 3.4e-7 * (300/Tn)^5 * exp(-3872/Tn)   # R25 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  HO2              -> H+(H2O)     + NO3            : 0.5e-9                                # R26 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  OH               -> H+(H2O)     + NO2            : 1.0e-10                               # R27 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  H                -> H+(H2O)     + NO             : 7.0e-12                               # R28 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  H2O        + M   -> NO+(H2O)2   + M              : 1.0e-27 * (308/Tn)^4.7                # R29 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  H2O        + M   -> NO+(H2O)3   + M              : 1.0e-27 * (308/Tn)^4.7                # R30 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)3  +  H2O              -> H+(H2O)3    + HONO           : 7.0e-11                               # R31 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)4   +  N2O5             -> H+(H2O)3(HNO3)  + HNO3       : 4.0e-12                               # R42 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)5   +  N2O5             -> H+(H2O)4(HNO3)  + HNO3       : 7.0e-12                               # R44 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)3(HNO3)     + H2O       -> H+(H2O)4    + HNO3           : 1.0e-9                                # R45 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)4(HNO3)     + H2O       -> H+(H2O)5    + HNO3           : 1.0e-9                                # R46 in Verronen et al. [2016] RNW")

# Verronen+ 2016: photochemistry in D region A2.2: recombination of positive ions with e- including RNW
reaction_rate_list.append(" NO+(N2)    +  e-    -> NO    + N2    : 1.4e-6 * (300/Tn)^0.4   # R5 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(CO2)   +  e-    -> NO    + CO2   : 1.5e-6                  # R6 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  e-    -> NO    + H2O   : 1.5e-6                  # R7 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  e-    -> NO    + 2H2O  : 2.0e-6                  # R8 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)3  +  e-    -> NO    + 3H2O  : 2.0e-6                  # R8 in Verronen et al. [2016] RNW")

# Verronen+ 2016: photochemistry in D region A2.5: negative ion reactions including RNW
reaction_rate_list.append(" O-         +  NO2          -> NO2-         + O            : 1.0e-9                             # R4  in Verronen et al. [2016] RNW")
reaction_rate_list.append(" O-         +  HNO3         -> NO3-         + OH           : 3.6e-9                             # R10 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" O2-        +  NO2          -> NO2-         + O2           : 7.0e-10                            # R14 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" O2-        +  HNO3         -> NO3-         + HO2          : 2.9e-9                             # R19 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" O3-        +  NO           -> NO3-         + O            : 1.05e-12 * (300/Tn)^2.15           # R23 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" O3-        +  NO2          -> NO3-         + O2           : 2.50e-11 * (300/Tn)^0.79           # R24 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" O3-        +  NO2          -> NO2-         + O3           : 7.5e-11 * (300/Tn)^0.79            # R25 in Verronen et al. [2016] RNW") # O3 in right side was O in original, but fixed as O3
reaction_rate_list.append(" O3-        +  NO           -> NO2-         + O2           : 1.05e-12 * (300/Tn)^2.15           # R26 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" CO3-       +  NO           -> NO2-         + CO2          : 1.3e-11 * (300/Tn)^1.64            # R38 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" CO3-       +  NO2          -> NO3-         + CO2          : 3.3e-11 * (300/Tn)^2.38            # R39 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" CO3-       +  HNO3         -> NO3-         + CO2    + OH  : 3.51e-10                           # R44 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO2-       +  H            -> OH-          + NO           : 3.0e-10                            # R51 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO2-       +  NO2          -> NO3-         + NO           : 2.0e-13                            # R52 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO2-       +  O3           -> NO3-         + O2           : 1.2e-10                            # R53 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO2-       +  H2O    + M   -> NO2-(H2O)    + M            : 1.6e-28                            # R57 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO2-       +  HNO3         -> NO3-         + HONO         : 1.6e-9                             # R58 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-       +  O            -> NO2-         + O2           : 0.5e-11                            # R59 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-       +  O3           -> NO2-         + 2O2          : 1.0e-13                            # R60 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-       +  H2O    + M   -> NO3-(H2O)    + M            : 1.6e-28                            # R62 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-       +  HNO3   + M   -> NO3-(HNO3)   + M            : 1.45e-26                           # R63 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" CO3-(H2O)  +  NO           -> NO2-         + H2O    + CO2 : 3.5e-12                            # R75 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" CO3-(H2O)  +  NO2          -> NO3-         + H2O    + CO2 : 4.0e-11                            # R76 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" CO3-(H2O)  +  NO2          -> NO3-(H2O)    + CO2          : 4.0e-11                            # R78 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" CO3-(H2O)  +  NO           -> NO2-(H2O)    + CO2          : 3.5e-12                            # R80 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO2-(H2O)  +  M            -> NO2-         + H2O    + M   : 5.7e-4 * (300/Tn) * exp(-7600/Tn)  # R82 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-(H2O)  +  H2O    + M   -> NO3-(H2O)2   + M            : 1.6e-28                            # R83 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-(H2O)  +  N2O5         -> NO3-(HNO3)   + HNO3         : 7.0e-10                            # R84 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-(H2O)  +  HNO3         -> NO3-(HNO3)   + H2O          : 1.6e-9                             # R85 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-(H2O)  +  M            -> NO3-         + H2O    + M   : 1.0e-3 * (300/Tn) * exp(-7300/Tn)  # R86 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-(H2O)2 +  M            -> NO3-(H2O)    + H2O    + M   : 1.5e-2 * (300/Tn) * exp(-7150/Tn)  # R87 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-(H2O)2 +  N2O5         -> NO3-(HNO3)   + HNO3   + H2O : 7.0e-10                            # R88 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-(HNO3) +  M            -> NO3-         + HNO3   + M   : 6.0e-3 * (300/Tn) * exp(-13130/Tn) # R90 in Verronen et al. [2016] RNW")

# Verronen+ 2016: photochemistry in D region A2.6: electron detachment from negative ions including RNW
reaction_rate_list.append(" O-     +  NO       -> NO2     + e-      : 3.1e-10 * (300/Tn)^0.83   # R2 in Verronen et al. [2016] RNW")

# Verronen+ 2016: photochemistry in D region A2.7: Photodetachment of electrons from negative ions including RNW
reaction_rate_list.append(" NO2-     -> NO2    + e-      : 8.0e-4    # R5 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO3-     -> NO3    + e-      : 5.2e-2    # R6 in Verronen et al. [2016] RNW")

# Verronen+ 2016: photochemistry in D region A2.9: ion-ion recombination two body reaction ions including RNW
reaction_rate_list.append(" H+(H2O)4   +  NO3-(HNO3)    -> 2HNO3  + 4H2O             : 6.0e-8 * (300/Tn)^0.5   # R1 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)4   +  NO3-          -> HNO3   + 4H2O             : 6.0e-8 * (300/Tn)^0.5   # R4 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)4   +  NO3-(H2O)     -> H      + 5H2O  + NO3      : 6.0e-8 * (300/Tn)^0.5   # R8 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)4   +  NO2-(H2O)     -> H      + 5H2O  + NO2      : 6.0e-8 * (300/Tn)^0.5   # R12 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)4   +  NO3-(H2O)2    -> H      + 6H2O  + NO3      : 6.0e-8 * (300/Tn)^0.5   # R15 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)4   +  NO2-          -> H      + 4H2O  + NO2      : 6.0e-8 * (300/Tn)^0.5   # R16 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)5   +  NO3-(HNO3)    -> 2HNO3  + 5H2O             : 6.0e-8 * (300/Tn)^0.5   # R17 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)5   +  NO3-          -> HNO3   + 5H2O             : 6.0e-8 * (300/Tn)^0.5   # R20 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)5   +  NO3-(H2O)     -> H      + 6H2O  + NO3      : 6.0e-8 * (300/Tn)^0.5   # R24 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)5   +  NO2-(H2O)     -> H      + 6H2O  + NO2      : 6.0e-8 * (300/Tn)^0.5   # R28 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)5   +  NO3-(H2O)2    -> H      + 7H2O  + NO3      : 6.0e-8 * (300/Tn)^0.5   # R31 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)5   +  NO2-          -> H      + 5H2O  + NO2      : 6.0e-8 * (300/Tn)^0.5   # R32 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)3   +  NO3-(HNO3)    -> 2HNO3  + 3H2O             : 6.0e-8 * (300/Tn)^0.5   # R33 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)3   +  NO3-          -> HNO3   + 3H2O             : 6.0e-8 * (300/Tn)^0.5   # R36 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)3   +  NO3-(H2O)     -> H      + 4H2O  + NO3      : 6.0e-8 * (300/Tn)^0.5   # R40 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)3   +  NO2-(H2O)     -> H      + 4H2O  + NO2      : 6.0e-8 * (300/Tn)^0.5   # R44 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)3   +  NO3-(H2O)2    -> H      + 5H2O  + NO3      : 6.0e-8 * (300/Tn)^0.5   # R47 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)3   +  NO2-          -> H      + 3H2O  + NO2      : 6.0e-8 * (300/Tn)^0.5   # R48 in Verronen et al. [2016] RNW")

reaction_rate_list.append(" NO+(H2O)   +  NO3-(HNO3)    -> NO  + H2O  + NO3  + HNO3  : 6.0e-8 * (300/Tn)^0.5   # R49 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  CO3-          -> NO  + H2O  + O    + CO2   : 6.0e-8 * (300/Tn)^0.5   # R50 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  NO3-          -> NO  + H2O  + NO3          : 6.0e-8 * (300/Tn)^0.5   # R52 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  HCO3-         -> NO  + H2O  + OH   + CO2   : 6.0e-8 * (300/Tn)^0.5   # R53 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  O2-           -> NO  + H2O  + O2           : 6.0e-8 * (300/Tn)^0.5   # R54 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  CO4-          -> NO  + H2O  + O2   + CO2   : 6.0e-8 * (300/Tn)^0.5   # R55 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  NO3-(H2O)     -> NO  + 2H2O + NO3          : 6.0e-8 * (300/Tn)^0.5   # R56 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  CO3-(H2O)2    -> NO  + 3H2O + O    + CO2   : 6.0e-8 * (300/Tn)^0.5   # R57 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  CO3-(H2O)     -> NO  + 2H2O + O    + CO2   : 6.0e-8 * (300/Tn)^0.5   # R59 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  NO2-(H2O)     -> NO  + 2H2O + NO2          : 6.0e-8 * (300/Tn)^0.5   # R60 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  NO3-(H2O)2    -> NO  + 3H2O + NO3          : 6.0e-8 * (300/Tn)^0.5   # R63 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)   +  NO2-          -> NO  + H2O  + NO2          : 6.0e-8 * (300/Tn)^0.5   # R64 in Verronen et al. [2016] RNW")

reaction_rate_list.append(" NO+(H2O)2  +  NO3-(HNO3)    -> NO  + 2H2O + NO3  + HNO3  : 6.0e-8 * (300/Tn)^0.5   # R65 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  CO3-          -> NO  + 2H2O + O    + CO2   : 6.0e-8 * (300/Tn)^0.5   # R66 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  NO3-          -> NO  + 2H2O + NO3          : 6.0e-8 * (300/Tn)^0.5   # R68 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  HCO3-         -> NO  + 2H2O + OH   + CO2   : 6.0e-8 * (300/Tn)^0.5   # R69 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  O2-           -> NO  + 2H2O + O2           : 6.0e-8 * (300/Tn)^0.5   # R70 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  CO4-          -> NO  + 2H2O + O2   + CO2   : 6.0e-8 * (300/Tn)^0.5   # R71 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  NO3-(H2O)     -> NO  + 3H2O + NO3          : 6.0e-8 * (300/Tn)^0.5   # R72 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  CO3-(H2O)2    -> NO  + 4H2O + O    + CO2   : 6.0e-8 * (300/Tn)^0.5   # R73 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  CO3-(H2O)     -> NO  + 3H2O + O    + CO2   : 6.0e-8 * (300/Tn)^0.5   # R75 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  NO2-(H2O)     -> NO  + 3H2O + NO2          : 6.0e-8 * (300/Tn)^0.5   # R76 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  NO3-(H2O)2    -> NO  + 4H2O + NO3          : 6.0e-8 * (300/Tn)^0.5   # R79 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+(H2O)2  +  NO2-          -> NO  + 2H2O + NO2          : 6.0e-8 * (300/Tn)^0.5   # R80 in Verronen et al. [2016] RNW")

reaction_rate_list.append(" NO+        +  NO3-(HNO3)    -> NO         + NO3  + HNO3  : 6.0e-8 * (300/Tn)^0.5   # R81 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  CO3-          -> NO         + O    + CO2   : 6.0e-8 * (300/Tn)^0.5   # R82 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  NO3-          -> NO         + NO3          : 6.0e-8 * (300/Tn)^0.5   # R84 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  HCO3-         -> NO         + OH   + CO2   : 6.0e-8 * (300/Tn)^0.5   # R85 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  O2-           -> NO         + O2           : 6.0e-8 * (300/Tn)^0.5   # R86 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  CO4-          -> NO         + O2   + CO2   : 6.0e-8 * (300/Tn)^0.5   # R87 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  NO3-(H2O)     -> NO  + H2O  + NO3          : 6.0e-8 * (300/Tn)^0.5   # R88 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  CO3-(H2O)2    -> NO  + 2H2O + O    + CO2   : 6.0e-8 * (300/Tn)^0.5   # R89 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  CO3-(H2O)     -> NO  + H2O  + O    + CO2   : 6.0e-8 * (300/Tn)^0.5   # R91 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  NO2-(H2O)     -> NO  + H2O  + NO2          : 6.0e-8 * (300/Tn)^0.5   # R92 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  NO3-(H2O)2    -> NO  + 2H2O + NO3          : 6.0e-8 * (300/Tn)^0.5   # R95 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" NO+        +  NO2-          -> NO         + NO2          : 6.0e-8 * (300/Tn)^0.5   # R96 in Verronen et al. [2016] RNW")

reaction_rate_list.append(" O2+        +  NO3-(HNO3)    -> O2         + NO3  + HNO3  : 6.0e-8 * (300/Tn)^0.5   # R97 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" O2+        +  NO3-          -> O2         + NO3          : 6.0e-8 * (300/Tn)^0.5   # R100 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" O2+        +  NO3-(H2O)     -> O2  + H2O  + NO3          : 6.0e-8 * (300/Tn)^0.5   # R104 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" O2+        +  NO2-(H2O)     -> O2  + H2O  + NO2          : 6.0e-8 * (300/Tn)^0.5   # R108 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" O2+        +  NO3-(H2O)2    -> O2  + 2H2O + NO3          : 6.0e-8 * (300/Tn)^0.5   # R111 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" O2+        +  NO2-          -> O2         + NO2          : 6.0e-8 * (300/Tn)^0.5   # R112 in Verronen et al. [2016] RNW")

# Verronen+ 2016: photochemistry in D region A2.10: ion-ion recombination three body reaction including RNW
reaction_rate_list.append(" H+(H2O)4   +  NO3-       + M   -> HNO3  + 4H2O   + M                : 1.25e-25 * (300/Tn)^4   # R2 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)5   +  NO3-       + M   -> HNO3  + 5H2O   + M                : 1.25e-25 * (300/Tn)^4   # R4 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)4   +  NO3-(HNO3) + M   -> 2HNO3 + 4H2O   + M                : 1.25e-25 * (300/Tn)^4   # R7 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)5   +  NO3-(HNO3) + M   -> 2HNO3 + 5H2O   + M                : 1.25e-25 * (300/Tn)^4   # R8 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)4   +  NO3-(H2O)  + M   -> H     + 5H2O   + NO3   + M        : 1.25e-25 * (300/Tn)^4   # R13 in Verronen et al. [2016] RNW")
reaction_rate_list.append(" H+(H2O)5   +  NO3-(H2O)  + M   -> H     + 6H2O   + NO3   + M        : 1.25e-25 * (300/Tn)^4   # R14 in Verronen et al. [2016] RNW")




# Chaffin+ 2017---------------------------------------------------------------------------------
reaction_rate_list.append(" O     + O   + M   -> O2   + M        : 5.4e-33 * (300/Tn)^3.25 # Chaffin et al. [2017] ")
reaction_rate_list.append(" O     + O2  + N2  -> O3   + N2       : 5.0e-35 * exp(724/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" O     + O2  + CO2 -> O3   + CO2      : 1.5e-33 * (300/Tn)^2.4 # Chaffin et al. [2017] ")
reaction_rate_list.append(" O     + O3        -> O2   + O2       : 8.0e-12 * exp(-2060/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" O     + CO  + M   -> CO2  + M        : 2.2e-33 * exp(-1780/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" O(1D) + O2        -> O    + O2       : 3.2e-11 * exp(70/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" O(1D) + O3        -> O2   + O2       : 1.2e-10 # Chaffin et al. [2017] ")
reaction_rate_list.append(" O(1D) + O3        -> O    + O  + O2  : 1.2e-10 # Chaffin et al. [2017] ")
reaction_rate_list.append(" O(1D) + H2        -> H    + OH       : 1.2e-10 # Chaffin et al. [2017] ")
reaction_rate_list.append(" O(1D) + CO2       -> O    + CO2      : 7.5e-11 * exp(115/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" O(1D) + H2O       -> OH   + OH       : 1.63e-10 * exp(60/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" H2    + O         -> OH   + H        : 6.34e-12 * exp(-4000/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" OH    + H2        -> H2O  + H        : 9.01e-13 * exp(-1526/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" H     + H   + CO2 -> H2   + CO2      : 1.6e-32 * (298/Tn)^2.27 # Chaffin et al. [2017] ")
reaction_rate_list.append(" H     + OH  + CO2 -> H2O  + CO2      : 1.292e-30 * (300/Tn)^2 # Chaffin et al. [2017] ")
reaction_rate_list.append(" H     + HO2       -> OH   + OH       : 7.2e-11 # Chaffin et al. [2017] ")
reaction_rate_list.append(" H     + HO2       -> H2O  + O(1D)    : 1.6e-12 # Chaffin et al. [2017] ")
reaction_rate_list.append(" H     + HO2       -> H2   + O2       : 3.45e-12 # Chaffin et al. [2017] ")
reaction_rate_list.append(" H     + H2O2      -> HO2  + H2       : 2.8e-12 * exp(-1890/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" H     + H2O2      -> H2O  + OH       : 1.7e-11 * exp(-1800/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" H     + O2  + M   -> HO2  + M        : k0 = 8.8e-32 * (300/Tn)^1.3 && \
                                                                   kinf = 7.5e-11 * (300/Tn)^(-0.2) # Chaffin et al. [2017] ")
reaction_rate_list.append(" H     + O3        -> OH   + O2       : 1.4e-10 * exp(-470/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" O     + OH        -> O2   + H        : 1.8e-11 * exp(180/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" O     + HO2       -> OH   + O2       : 3.0e-11 * exp(200/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" O     + H2O2      -> OH   + HO2      : 1.4e-12 * exp(-2000/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" OH    + OH        -> H2O  + O        : 1.8e-12 # Chaffin et al. [2017] ")
reaction_rate_list.append(" OH    + OH  + M   -> H2O2 + M        : k0 = 8.97e-31 * (300/Tn)^1 && \
                                                                   kinf = 2.6e-11 # Chaffin et al. [2017] ")
reaction_rate_list.append(" OH    + O3        -> HO2  + O2       : 1.7e-12 * exp(-940/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" OH    + HO2       -> H2O  + O2       : 4.8e-11 * exp(250/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" OH    + H2O2      -> H2O  + HO2      : 1.8e-12 # Chaffin et al. [2017] ")
reaction_rate_list.append(" HO2   + O3        -> OH   + O2 + O2  : 1.0e-14 * exp(-490/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" HO2   + HO2       -> H2O2 + O2       : 3.0e-13 * exp(460/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" HO2   + HO2 + M   -> H2O2 + O2 + M   : 4.2e-33 * exp(920/Tn) # Chaffin et al. [2017] ")
reaction_rate_list.append(" CO    + OH  + M   -> CO2  + H  + M   : k0M = 1.5e-13 * (300/Tn)^(-0.6) && \
                                                                   kinfM = 2.1e9 * (300/Tn)^(-6.1) # Chaffin et al. [2017] ")
reaction_rate_list.append(" CO    + OH  + M   -> HOCO + M        : k0 = 5.9e-33 * (300/Tn)^1.4 && \
                                                                   kinf = 1.1e-12 * (300/Tn)^(-1.3) # Chaffin et al. [2017] ")
reaction_rate_list.append(" HOCO  + O2        -> HO2  + CO2      : 2.0e-12 # Chaffin et al. [2017] ")
reaction_rate_list.append(" CO2+  + H2        -> CO2  + H  + H   : 8.7e-10 # Chaffin et al. [2017] ") # it's only used to reproduce Chaffin+2017, otherwise remove it!
#reaction_rate_list.append(" CO    + OH        -> CO2  + H        : 1.5e-13 # Chaffin et al. [2017] ")

#Yoshida et al. [2023]
reaction_rate_list.append(" ^13CO2   + hv -> ^13CO     + O           : Photodissociation # Yoshida et al. [2023] ")
reaction_rate_list.append(" ^13CO2   + hv -> ^13CO     + O(1D)       : Photodissociation # Yoshida et al. [2023] ")
reaction_rate_list.append(" O     + ^13CO + M -> ^13CO2  + M     : 1.0074 * 2.2e-33 * exp(-1780/Tn) # Yoshida et al. [2023] ")
reaction_rate_list.append(" ^13CO + OH + M    -> ^13CO2  + H  + M: k0M = 0.98911 * 1.5e-13 * (300/Tn)^(-0.6) && \
                                                                   kinfM = 0.98911 * 2.1e9 * (300/Tn)^(-6.1) # Yoshida et al. [2023] ")
reaction_rate_list.append(" ^13CO + OH + M    -> HO^13CO + M     : k0 = 0.98911 * 5.9e-33 * (300/Tn)^1.4 && \
                                                                   kinf = 0.98911 * 1.1e-12 * (300/Tn)^(-1.3) # Yoshida et al. [2023] ")
reaction_rate_list.append(" HO^13CO  + O2     -> HO2  + ^13CO2   : 0.98911*2.0e-12 # Yoshida et al. [2023] ")


## NOx reactions JPL document ---------------------------------------------------------------------------------
#reaction_rate_list.append(" O     + NO    + M   -> NO2      + M        : k0 = 9.0e-32 * (Tn/300)^(-1.5) && \
#                                                                         kinf = 3.0e-11 ")
#reaction_rate_list.append(" O     + NO2         -> NO       + O2       : 5.1e-12 * exp(210/Tn) ")
#reaction_rate_list.append(" O     + NO2   + M   -> NO3      + M        : k0 = 2.5e-31 * (Tn/300)^(-1.8) && \
#                                                                         kinf = 2.2e-11 * (Tn/300)^(-0.7) ")
#reaction_rate_list.append(" O     + NO3         -> O2       + NO2      : 1.0e-11 ")
#reaction_rate_list.append(" O     + N2O5        -> products            : 0.0 ")
#reaction_rate_list.append(" O     + HNO3        -> OH       + NO3      : 0.0 ")
#reaction_rate_list.append(" O     + HO2NO2      -> products            : 7.8e-11 * exp(-3400/Tn) ")
#reaction_rate_list.append(" H     + NO2         -> OH       + NO       : 4.0e-10 * exp(-340/Tn) ")
#reaction_rate_list.append(" OH    + NO    + M   -> HONO     + M        : k0 = 7.0e-31 * (Tn/300)^(-2.6) && \
#                                                                         kinf = 3.6e-11 * (Tn/300)^(-0.1) ")
#reaction_rate_list.append(" OH    + NO2   + M   -> HNO3     + M        : k0 = 1.8e-30 * (Tn/300)^(-3.0) && \
#                                                                         kinf = 2.8e-11 ") # HNO3 = HONO2
#reaction_rate_list.append(" OH    + NO3         -> O2       + NO2      : 2.2e-11 ") # @ 298 K
#reaction_rate_list.append(" OH    + HONO        -> H2O      + NO2      : 1.8e-11 * exp(-390/Tn) ")
#reaction_rate_list.append(" OH    + HNO3        -> H2O      + NO3      : k0LH = 2.4e-14 * exp(460/Tn) && \
#                                                                         k2LH = 2.7e-17 * exp(2199/Tn) && \
#                                                                         k3LH = 6.5e-34 * exp(1335/Tn) ") #see note, LH: Lindemann-Hinshelwood expression
#reaction_rate_list.append(" OH    + HO2NO2      -> H2O      + NO2 + O2 : 1.3e-12 * exp(380/Tn) ") # ref: Nair et al., 1994
#reaction_rate_list.append(" OH    + NH3         -> H2O      + NH2      : 1.7e-12 * exp(-710/Tn) ")
#reaction_rate_list.append(" HO2   + NO          -> NO2      + OH       : 3.3e-12 * exp(270/Tn) ")
#reaction_rate_list.append(" NO2*  + H2O         -> OH       + HONO     : 0.0 ")# see note
#reaction_rate_list.append(" HO2   + NO2   + M   -> HO2NO2   + M        : k0 = 1.9e-31 * (Tn/300)^(-3.4) && \
#                                                                         kinf = 4.0e-12 * (Tn/300)^(-0.3) ")
#reaction_rate_list.append(" HO2   + NO2         -> HONO     + O2       : 0.0 ")#see note
#reaction_rate_list.append(" HO2   + NO3         -> products            : 3.5e-12 ") # @ 298 K
#reaction_rate_list.append(" HO2   + NH2         -> products            : 3.4e-11 ")
#reaction_rate_list.append(" N     + O2          -> NO       + O        : 1.5e-11 * exp(-3600/Tn) ")
#reaction_rate_list.append(" N     + O3          -> NO       + O2       : 0.0 ")
#reaction_rate_list.append(" N     + NO          -> N2       + O        : 2.1e-11 * exp(100/Tn) ")
#reaction_rate_list.append(" N     + NO2         -> N2O      + O        : 5.8e-12 * exp(220/Tn) ")
#reaction_rate_list.append(" NO    + O3          -> NO2      + O2       : 3.0e-12 * exp(-1500/Tn) ")
#reaction_rate_list.append(" NO    + NO3         -> 2NO2                : 1.5e-11 * exp(170/Tn) ")
#reaction_rate_list.append(" NO2   + O3          -> NO3      + O2       : 1.2e-13 * exp(-2450/Tn) ")
#reaction_rate_list.append(" NO2   + NO3         -> NO       + NO2 + O2 : 4.5e-14 * exp (-1260/Tn) ")#see note
#reaction_rate_list.append(" NO2   + NO3   + M   -> N2O5     + M        : k0 = 2.4e-30 * (Tn/300)^(-3.0) && \
#                                                                         kinf = 1.6e-12 * (Tn/300)^(0.1) ")
#reaction_rate_list.append(" NO3   + NO3         -> 2NO2     + O2       : 8.5e-13 * exp(-2450/Tn) ")
#reaction_rate_list.append(" NH2   + O2          -> products            : 0.0 ")
#reaction_rate_list.append(" NH2   + O3          -> products            : 4.3e-12 * exp(-930/Tn) ")
#reaction_rate_list.append(" NH2   + NO          -> products            : 4.0e-12 * exp(450/Tn) ")
#reaction_rate_list.append(" NH2   + NO2         -> products            : 2.1e-12 * exp(650/Tn) ")
#reaction_rate_list.append(" NH    + NO          -> products            : 4.9e-11 ")
#reaction_rate_list.append(" NH    + NO2         -> products            : 3.5e-13 * exp(1140/Tn) ")
#reaction_rate_list.append(" O3    + HONO        -> O2       + HNO3     : 0.0 ")
#reaction_rate_list.append(" N2O5  + H2O         -> 2HNO3               : 0.0 ")
#
## Production and Loss of N2O
## JPL
#reaction_rate_list.append(" O(1D) + N2     + M  -> N2O      + M        : 2.8e-36 * (Tn/300)^(-0.9) ")
## Adams et al., 2021, Airapetian et al., 2016
#reaction_rate_list.append(" N2O   + CH          -> NO       + HCN      : 1.5e-11 * exp(257/Tn) ")
#reaction_rate_list.append(" N2O   + N(2D)       -> N2       + NO       : 1.5e-11 * exp(-570/Tn) ")
## JPL
#reaction_rate_list.append(" N2O   + O(1D)       -> N2       + O2       : 0.39 * 1.9e-10 * exp(20/Tn) ")
#reaction_rate_list.append(" N2O   + O(1D)       -> NO       + NO       : 0.61 * 1.9e-10 * exp(20/Tn) ")
#
## Loss of N : Adams et al., 2021
#reaction_rate_list.append(" N     + N      + M  -> N2       + M        : 8.27e-34 * exp(490/Tn) ")
#reaction_rate_list.append(" N     + O      + M  -> NO       + M        : 5.46e-33 * exp(155/Tn) ")
#
#
### Krasnopolsky, 2013
##reaction_rate_list.append(" N     + O           -> NO                  : 3.3e-16 * Tn^(-0.5) # Krasnopolsky [2013] ")
##reaction_rate_list.append(" N     + O      + M  -> NO       + M        : 3.5e-31 * Tn^(-0.5) # Krasnopolsky [2013] ")
##reaction_rate_list.append(" N     + NO          -> N2       + O        : 2.1e-11 * exp(100/Tn) # Krasnopolsky [2013] ")
##reaction_rate_list.append(" NO    + O      + M  -> NO2      + M        : 1.0e-27 * Tn^(-1.5) # Krasnopolsky [2013] ")
##reaction_rate_list.append(" NO    + O3          -> NO2      + O2       : 3.0e-12 * exp(-1500/Tn) # Krasnopolsky [2013] ")
##reaction_rate_list.append(" NO    + HO2         -> NO2      + OH       : 3.5e-12 * exp(250/Tn) # Krasnopolsky [2013] ")
##reaction_rate_list.append(" NO2   + O           -> NO       + O2       : 5.1e-12 * exp(210/Tn) # Krasnopolsky [2013] ")
##reaction_rate_list.append(" NO2   + O      + M  -> NO3      + M        : 1.7e-26 * Tn^(-1.8) # Krasnopolsky [2013] ")
##reaction_rate_list.append(" NO3   + NO          -> 2NO2                : 1.5e-11 * exp(170/Tn) # Krasnopolsky [2013] ")
##reaction_rate_list.append(" NO3   + O           -> NO2      + O2       : 1.0e-11  # Krasnopolsky [2013] ")
##reaction_rate_list.append(" N     + O2          -> NO       + O        : 2.0e-18 Tn^2.15 * exp(-25600/Tn) # Krasnopolsky [2013] ")
#
## Krasnopolsky, 2006
#reaction_rate_list.append(" NO    + NO          -> N2O      + O        : 6e-12 * exp(-32900/Tn) # Krasnopolsky [2006]")


# Nair et al., 1994: NOx reactions
reaction_rate_list.append(" N      + O2          -> NO       + O        : 1.5e-11 * exp(-3600/Tn)  # R57 in Nair et al. [1994] ")
reaction_rate_list.append(" N      + O3          -> NO       + O2       : 1.0e-16  # R58 in Nair et al. [1994] ")
reaction_rate_list.append(" N      + OH          -> NO       + H        : 3.8e-11 * exp(85/Tn)  # R59 in Nair et al. [1994] ")
reaction_rate_list.append(" N      + HO2         -> NO       + OH       : 2.2e-11  # R60 in Nair et al. [1994] ")
reaction_rate_list.append(" N      + NO          -> N2       + O        : 3.4e-11  # R61 in Nair et al. [1994] ")
reaction_rate_list.append(" N      + NO2         -> N2O      + O        : 3.0e-12  # R62 in Nair et al. [1994] ")
reaction_rate_list.append(" N(2D)  + O           -> N        + O        : 6.9e-13  # R63 in Nair et al. [1994] ")
reaction_rate_list.append(" N(2D)  + CO2         -> NO       + CO       : 3.5e-13  # R64 in Nair et al. [1994] ")
reaction_rate_list.append(" N(2D)  + N2          -> N        + N2       : 1.7e-14  # R65 in Nair et al. [1994] ")
reaction_rate_list.append(" N(2D)  + NO          -> N2       + O        : 6.9e-11  # R66 in Nair et al. [1994] ")
reaction_rate_list.append(" O      + NO    + M   -> NO2      + M        : k0   = 1.2e-27 * Tn^(-1.5) && \
                                                                          kinf = 3.0e-11  # R67 in Nair et al. [1994] ")
reaction_rate_list.append(" O      + NO2         -> NO       + O2       : 6.5e-12 * exp(120/Tn)  # R68 in Nair et al. [1994] ")
reaction_rate_list.append(" O      + NO2   + M   -> NO3      + M        : k0   = 2.0e-26 * Tn^(-2.0) && \
                                                                          kinf = 2.2e-11  # R69 in Nair et al. [1994] ")
reaction_rate_list.append(" O      + NO3         -> O2       + NO2      : 1.0e-11  # R70 in Nair et al. [1994] ")
reaction_rate_list.append(" O      + HO2NO2      -> OH       + NO2 + O2 : 7.8e-11 * exp(-3400/Tn)  # R71 in Nair et al. [1994] ")
reaction_rate_list.append(" O(1D)  + N2          -> O        + N2       : 1.8e-11 * exp(110/Tn)  # R72 in Nair et al. [1994] ")
reaction_rate_list.append(" O(1D)  + N2    + M   -> N2O      + M        : 2.8e-35 * Tn^(-0.6)  # R73 in Nair et al. [1994] ")
reaction_rate_list.append(" O(1D)  + N2O         -> 2NO                 : 6.7e-11  # R74 in Nair et al. [1994] ")
reaction_rate_list.append(" O(1D)  + N2O         -> N2       + O2       : 4.9e-11  # R75 in Nair et al. [1994] ")
reaction_rate_list.append(" NO     + O3          -> NO2      + O2       : 2.0e-12 * exp(-1400/Tn)  # R76 in Nair et al. [1994] ")
reaction_rate_list.append(" NO     + HO2         -> NO2      + OH       : 3.7e-12 * exp(240/Tn)  # R77 in Nair et al. [1994] ")
reaction_rate_list.append(" NO     + NO3         -> 2NO2                : 1.7e-11 * exp(150/Tn)  # R78 in Nair et al. [1994] ")
reaction_rate_list.append(" H      + NO2         -> OH       + NO       : 2.2e-10 * exp(-182/Tn)  # R79 in Nair et al. [1994] ")
reaction_rate_list.append(" H      + NO3         -> OH       + NO2      : 1.1e-10  # R80 in Nair et al. [1994] ")
reaction_rate_list.append(" OH     + NO    + M   -> HONO     + M        : k0   = 4.8e-24 * Tn^(-2.6) && \
                                                                          kinf = 2.6e-10 * Tn^(-0.5)  # R81 in Nair et al. [1994] ")
reaction_rate_list.append(" OH     + NO2   + M   -> HNO3     + M        : k0   = 5.5e-22 * Tn^(-3.2) && \
                                                                          kinf = 4.0e-8 * Tn^(-1.3)  # R82 in Nair et al. [1994] ")
reaction_rate_list.append(" OH     + NO3         -> HO2      + NO2      : 2.3e-11  # R83 in Nair et al. [1994] ")
reaction_rate_list.append(" OH     + HONO        -> H2O      + NO2      : 1.8e-11 * exp(-390/Tn)  # R84 in Nair et al. [1994] ")
reaction_rate_list.append(" OH     + HNO3        -> H2O      + NO3      : 7.2e-15 * exp(785/Tn)  # R85 in Nair et al. [1994] ")
reaction_rate_list.append(" OH     + HO2NO2      -> H2O      + NO2 + O2 : 1.3e-12 * exp(380/Tn)  # R86 in Nair et al. [1994] ")
reaction_rate_list.append(" HO2    + NO2   + M   -> HO2NO2   + M        : k0   = 3.8e-23 * Tn^(-3.2) && \
                                                                          kinf = 1.4e-8 * Tn^(-1.4)  # R87 in Nair et al. [1994] ")
reaction_rate_list.append(" HO2    + NO3         -> O2       + HNO3     : 9.2e-13  # R88 in Nair et al. [1994] ")
reaction_rate_list.append(" NO2    + O3          -> NO3      + O2       : 1.2e-13 * exp(-2450/Tn)  # R89 in Nair et al. [1994] ")
reaction_rate_list.append(" NO2    + NO3   + M   -> N2O5     + M        : k0   = 2.5e-19 * Tn^(-4.3) && \
                                                                          kinf = 2.6e-11 * Tn^(-0.5)  # R90 in Nair et al. [1994] ")
reaction_rate_list.append(" NO2    + NO3         -> NO       + NO2 + O2 : 8.2e-14 * exp(-1480/Tn)  # R91 in Nair et al. [1994] ")

# N + CO2 -> NO + CO   : Rawlins and Kaufman [1976]
reaction_rate_list.append(" N     + CO2          -> NO    + CO      : 1.0e-19 @ upper limit at 300 K by Rawlins and Kaufman [1976] ")


# Pinto et al. [1980] Formaldehyde H2CO reactions
reaction_rate_list.append(" H     + CO     + M   -> HCO   + M      : 2.0e-33 * exp(-850/Tn) # R7 in Pinto et al. [1980] ")
reaction_rate_list.append(" H     + HCO          -> H2    + CO     : 3.0e-10 # R8 in Pinto et al. [1980] ")
reaction_rate_list.append(" HCO   + HCO          -> H2CO  + CO     : 6.3e-11 # R9 in Pinto et al. [1980] ")
reaction_rate_list.append(" OH    + HCO          -> H2O   + CO     : 5.0e-11 # R10 in Pinto et al. [1980] ")
reaction_rate_list.append(" O     + HCO          -> H     + CO2    : 1.0e-10 # R11 in Pinto et al. [1980] ")
reaction_rate_list.append(" O     + HCO          -> OH    + CO     : 1.0e-10 # R12 in Pinto et al. [1980] ")
reaction_rate_list.append(" O2    + HCO          -> HO2   + CO     : 5.0e-12 # R13 in Pinto et al. [1980] ")
reaction_rate_list.append(" HO2   + HCO          -> H2O2  + CO     : 1.0e-11 # R14 in Pinto et al. [1980] ")
reaction_rate_list.append(" H     + H2CO         -> H2    + HCO    : 2.8e-11 * exp(-1540/Tn) # R15 in Pinto et al. [1980] ")
reaction_rate_list.append(" OH    + H2CO         -> H2O   + HCO    : 1.7e-11 * exp(-100/Tn) # R16 in Pinto et al. [1980] ")

reaction_rate_list.append(" H2CO                 -> Rainout        : Rainout ")

# --------------------------------------------------------------------------------------------------------------------------
#
#                                                       Jupiter
#
# --------------------------------------------------------------------------------------------------------------------------
Planet_list.append(['Jupiter',len(reaction_rate_list)])
# Photoionization
reaction_rate_list.append(" H2   + e-* -> H2+   + e- + e-*      : Impact ionization ")
reaction_rate_list.append(" H    + hv  -> H+    + e-            : Photoionization ")
reaction_rate_list.append(" H2   + hv  -> H2+   + e-            : Photoionization ")
reaction_rate_list.append(" H2   + hv  -> H+    + e- + H        : Photoionization ")
reaction_rate_list.append(" He   + hv  -> He+   + e-            : Photoionization ")
reaction_rate_list.append(" CH4  + hv  -> CH4+  + e-            : Photoionization ")
reaction_rate_list.append(" CH4  + hv  -> CH3+  + e- + H        : Photoionization ")
reaction_rate_list.append(" CH4  + hv  -> H2+   + e- + products : Photoionization ")
reaction_rate_list.append(" C2H2 + hv  -> C2H2+ + e-            : Photoionization ")
reaction_rate_list.append(" C2H4 + hv  -> C2H4+ + e-            : Photoionization ")
reaction_rate_list.append(" C2H4 + hv  -> C2H3+ + e- + H        : Photoionization ")
reaction_rate_list.append(" C2H4 + hv  -> C2H2+ + e- + products : Photoionization ")
reaction_rate_list.append(" C2H4 + hv  -> C2H+  + e- + products : Photoionization ")
reaction_rate_list.append(" C2H6 + hv  -> C2H6+ + e-            : Photoionization ")
reaction_rate_list.append(" C2H6 + hv  -> C2H5+ + e- + H        : Photoionization ")
reaction_rate_list.append(" C2H6 + hv  -> C2H4+ + e- + products : Photoionization ")
reaction_rate_list.append(" C2H6 + hv  -> C2H3+ + e- + products : Photoionization ")
reaction_rate_list.append(" C2H6 + hv  -> C2H2+ + e- + products : Photoionization ")

# Kim and Fox 1994 : Recombination
reaction_rate_list.append(" H+    + e- -> H               : 4.0e-12* ( Te/250.0)^(-0.7) @ Bates and Dalgarno [1962] # RC1 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + e- -> He              : 4.0e-12* ( Te/250.0)^(-0.7) @ Bates and Dalgarno [1962] # RC2 in Kim and Fox [1994] ")
reaction_rate_list.append(" HeH+  + e- -> H    +  He      : 1.0e-8 * ( Te/300.0)^(-0.6) @ Yousif and Mitchell [1989] # RC3 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + e- -> H    +  H       : 2.3e-7 * ( Te/300.0)^(-0.4) @ Auerbach et al. [1977] # RC4 in Kim and Fox [1994] ")
reaction_rate_list.append(" H3+   + e- -> H2   +  H       : 4.4e-8 * ( Te/300.0)^(-0.5) @ Canosa et al. [1992] # RC5 in Kim and Fox [1994] ")
reaction_rate_list.append(" H3+   + e- -> H    +  H +  H  : 5.6e-8 * ( Te/300.0)^(-0.5) @ Mitchell et al. [1983] # RC6 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + e- -> C               : 4.0e-12* ( Te/250.0)^(-0.7) @ Bates and Dalgarno [1962] # RC7 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH+   + e- -> C    +  H       : 1.5e-7 * ( Te/300.0)^(-0.42) @ Mitchell and McGowan [1978] # RC8 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH2+  + e- -> CH   +  H       : 2.5e-7 * ( Te/300.0)^(-0.5) @ Mul et al. [1981] # RC9 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH3+  + e- -> CH2  +  H       : 3.5e-7 * ( Te/300.0)^(-0.5) @ Mul et al. [1981] # RC10 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH4+  + e- -> CH3  +  H       : 3.5e-7 * ( Te/300.0)^(-0.5) @ Mul et al. [1981] # RC11 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH4+  + e- -> CH2  +  H +  H  : 3.5e-7 * ( Te/300.0)^(-0.5) @ Mul et al. [1981] # RC12 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH5+  + e- -> CH2  +  H +  H2 : 8.8e-7 * ( Te/300.0)^(-0.5) @ Adams et al. [1984] # RC13 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH5+  + e- -> CH3  +  H +  H  : 2.2e-7 * ( Te/300.0)^(-0.5) @ Adams et al. [1991] # RC14 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2+   + e- -> C    +  C       : 3.0e-7 * ( Te/300.0)^(-0.5) @ Mul and McGowan [1980] # RC15 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H+  + e- -> C2   +  H       : 2.7e-7 * ( Te/300.0)^(-0.5) @ Mul and McGowan [1980] # RC16 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H+  + e- -> CH   +  C       : 2.7e-7 * ( Te/300.0)^(-0.5) @ Mul and McGowan [1980] # RC17 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + e- -> C2H  +  H       : 2.7e-7 * ( Te/300.0)^(-0.5) @ Mul and McGowan [1980] # RC18 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + e- -> CH   +  CH      : 2.7e-7 * ( Te/300.0)^(-0.5) @ Mul and McGowan [1980] # RC19 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H3+ + e- -> C2H2 +  H       : 4.5e-7 * ( Te/300.0)^(-0.5) @ Mul and McGowan [1980] # RC20 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H3+ + e- -> CH2  +  CH      : 4.5e-7 * ( Te/300.0)^(-0.5) @ Mul and McGowan [1980] # RC21 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H4+ + e- -> C2H3 +  H       : 3.0e-7 * ( Te/300.0)^(-0.5) @ estimate # RC22 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H4+ + e- -> CH2  +  CH2     : 3.0e-7 * ( Te/300.0)^(-0.5) @ estimate # RC23 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H5+ + e- -> C2H4 +  H       : 7.4e-7 * ( Te/300.0)^(-0.5) @ Adams and Smith [1988] # RC24 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H5+ + e- -> CH3  +  CH2     : 7.4e-7 * ( Te/300.0)^(-0.5) @ Adams and Smith [1988] # RC25 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H6+ + e- -> C2H5 +  H       : 3.0e-7 * ( Te/300.0)^(-0.5) @ estimate # RC26 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H6+ + e- -> CH3  +  CH3     : 3.0e-7 * ( Te/300.0)^(-0.5) @ estimate # RC27 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H7+ + e- -> C2H6 +  H       : 3.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC28 in Kim and Fox [1994] ")
reaction_rate_list.append(" C3H+  + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC29 in Kim and Fox [1994] ")
reaction_rate_list.append(" C3H2+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC29 in Kim and Fox [1994] ")
reaction_rate_list.append(" C3H3+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC29 in Kim and Fox [1994] ")
reaction_rate_list.append(" C3H4+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC29 in Kim and Fox [1994] ")
reaction_rate_list.append(" C3H5+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC29 in Kim and Fox [1994] ")
reaction_rate_list.append(" C3H6+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC29 in Kim and Fox [1994] ")
reaction_rate_list.append(" C3H7+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC29 in Kim and Fox [1994] ")
reaction_rate_list.append(" C3H8+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC29 in Kim and Fox [1994] ")
reaction_rate_list.append(" C3H9+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC29 in Kim and Fox [1994] ")
reaction_rate_list.append(" C4H+  + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC30 in Kim and Fox [1994] ")
reaction_rate_list.append(" C4H2+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC30 in Kim and Fox [1994] ")
reaction_rate_list.append(" C4H3+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC30 in Kim and Fox [1994] ")
reaction_rate_list.append(" C4H5+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC30 in Kim and Fox [1994] ")
reaction_rate_list.append(" C4H7+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC30 in Kim and Fox [1994] ")
reaction_rate_list.append(" C4H9+ + e- -> products        : 7.5e-7 * ( Te/300.0)^(-0.5) @ estimate # RC30 in Kim and Fox [1994] ")

# Kim and Fox 1994 : neutral - ions
reaction_rate_list.append(" H2+   + H2       -> H3+   + H             : 2.00e-9 @ Clow and Futrell [1972] # R1 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + H        -> H+    + H2            : 6.40e-10 @ Karpas et al. [1979] # R2 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + He       -> HeH+  + H             : 1.40e-10 @ Smith et al. [1976] # R3 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + CH4      -> CH5+  + H             : 1.10e-10 @ Kim and Huntress [1975b] # R4 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + CH4      -> CH4+  + H2            : 1.41e-9 @ Kim and Huntress [1975b] # R5 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + CH4      -> CH3+  + H   + H2      : 2.28e-9 @ Kim and Huntress [1975b] # R6 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + C2H2     -> C2H3+ + H             : 4.77e-10 @ Kim and Huntress [1975b] # R7 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + C2H2     -> C2H2+ + H2            : 4.82e-9 @ Kim and Huntress [1975b] # R8 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + C2H6     -> C2H6+ + H2            : 2.94e-10 @ Kim and Huntress [1975b] # R9 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + C2H6     -> C2H5+ + H   + H2      : 1.37e-9 @ Kim and Huntress [1975b] # R10 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + C2H6     -> C2H4+ + H2  + H2      : 2.35e-9 @ Kim and Huntress [1975b] # R11 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + C2H6     -> C2H3+ + H   + H2 + H2 : 6.86e-10 @ Kim and Huntress [1975b] # R12 in Kim and Fox [1994] ")
reaction_rate_list.append(" H2+   + C2H6     -> C2H2+ + H2  + H2 + H2 : 1.96e-10 @ Kim and Huntress [1975b] # R13 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + H2       -> H2+   + He            : 9.35e-15 @ Schauer et al. [1989] # R14 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + H2       -> H+    + H   + He      : 4.57e-14 @ Schauer et al. [1989] # R15 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + H2(v>=2) -> H+    + H   + He      : 1.0e-9  @ Jones et al. [1986] # R16 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + CH4      -> H+    + CH3 + He      : 4.76e-10 @ Mauclaire et al. [1978] # R17 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + CH4      -> CH+   + H2  + H  + He : 2.38e-10 @ Mauclaire et al. [1978] # R18 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + CH4      -> CH2+  + H2  + He      : 8.50e-10 @ Mauclaire et al. [1978] # R19 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + CH4      -> CH3+  + H   + He      : 8.50e-11 @ Mauclaire et al. [1978] # R20 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + CH4      -> CH4+  + He            : 5.10e-11 @ Mauclaire et al. [1978] # R21 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H2     -> C2H2+ + He            : 2.45e-10 @ Kim and Huntress [1975a] # R22 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H2     -> C2H+  + H   + He      : 8.75e-10 @ Kim and Huntress [1975a] # R23 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H2     -> C2+   + H2  + He      : 1.61e-9 @ Kim and Huntress [1975a] # R24 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H2     -> CH+   + CH  + He      : 7.70e-10 @ Kim and Huntress [1975a] # R25 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H4     -> C2H4+ + He            : 2.38e-10 @ Kim and Huntress [1975a] # R26 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H4     -> C2H3+ + H   + He      : 1.70e-10 @ Kim and Huntress [1975a] # R27 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H4     -> C2H2+ + H2  + He      : 2.18e-9 @ Kim and Huntress [1975a] # R28 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H4     -> C2H+  + H   + H2 + He : 4.42e-10 @ Kim and Huntress [1975a] # R29 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H4     -> CH2+  + CH2 + He      : 4.08e-10 @ Kim and Huntress [1975a] # R30 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H6     -> C2H4+ + H2  + He      : 4.20e-10 @ Kim and Huntress [1975a] # R31 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H6     -> C2H3+ + H   + H2 + He : 1.74e-9 @ Kim and Huntress [1975a] # R32 in Kim and Fox [1994] ")
reaction_rate_list.append(" He+   + C2H6     -> C2H2+ + H2  + H2 + He : 8.40e-10 @ Kim and Huntress [1975a] # R33 in Kim and Fox [1994] ")
reaction_rate_list.append(" H+    + CH4      -> CH3+  + H2            : 3.69e-9 @ Dheandhanoo et al. [1984] # R34 in Kim and Fox [1994] ")
reaction_rate_list.append(" H+    + CH4      -> CH4+  + H             : 0.81e-9 @ Dheandhanoo et al. [1984] # R35 in Kim and Fox [1994] ")
reaction_rate_list.append(" H+    + C2H6     -> C2H3+ + H2 +  H2      : 1.30e-9 @ Mackay et al. [1981] # R36 in Kim and Fox [1994] ")
reaction_rate_list.append(" H+    + C2H6     -> C2H4+ + H  +  H2      : 1.30e-9 @ Mackay et al. [1981] # R37 in Kim and Fox [1994] ")
reaction_rate_list.append(" H+    + C2H6     -> C2H5+ + H2            : 1.30e-9 @ Mackay et al. [1981] # R38 in Kim and Fox [1994] ")
reaction_rate_list.append(" H+    + H2(v>=4) -> H2+   + H             : 2.0e-9 # Kim and Fox [1994]")
reaction_rate_list.append(" HeH+  + H2       -> H3+   + He            : 1.50e-9 @ Bohme et al. [1980] # R39 in Kim and Fox [1994] ")
reaction_rate_list.append(" HeH+  + H        -> H2+   + He            : 9.10e-10 @ Karpas et al. [1979] # R40 in Kim and Fox [1994] ")
reaction_rate_list.append(" HeH+  + C2H4     -> C2H4+ + H  +  He      : 7.00e-10 @ Smith and Futrell [1976] # R41 in Kim and Fox [1994] ")
reaction_rate_list.append(" HeH+  + C2H4     -> C2H3+ + H2 +  He      : 2.10e-9 @ Smith and Futrell [1976] # R42 in Kim and Fox [1994] ")
reaction_rate_list.append(" HeH+  + C2H6     -> C2H3+ + H2 +  H2 + He : 1.05e-9 @ Mackay et al. [1981] # R43 in Kim and Fox [1994] ")
reaction_rate_list.append(" HeH+  + C2H6     -> C2H5+ + H2 +  He      : 1.05e-9 @ Mackay et al. [1981] # R44 in Kim and Fox [1994] ")
reaction_rate_list.append(" H3+   + CH4      -> CH5+  + H2            : 2.40e-9 @ Bohme et al. [1980] # R45 in Kim and Fox [1994] ")
reaction_rate_list.append(" H3+   + C2H2     -> C2H3+ + H2            : 2.90e-9 @ Mackay et al. [1977] # R46 in Kim and Fox [1994] ")
reaction_rate_list.append(" H3+   + C2H4     -> C2H5+ + H2            : 0.69e-9 @ Rakshit [1982] # R47 in Kim and Fox [1994] ")
reaction_rate_list.append(" H3+   + C2H4     -> C2H3+ + H2  +  H2     : 1.61e-9 @ Rakshit [1982] # R48 in Kim and Fox [1994] ")
reaction_rate_list.append(" H3+   + C2H6     -> C2H5+ + H2  +  H2     : 2.40e-9 @ Mackay et al. [1981] # R49 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + CH4      -> C2H2+ + H2            : 3.30e-10 @ Bohme et al. [1982] # R50 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + CH4      -> C2H3+ + H             : 9.80e-10 @ Bohme et al. [1982] # R51 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + C2H2     -> C3H+  + H             : 2.80e-9 @ Anicich et al. [1986] # R52 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + C2H4     -> C2H3+ + CH            : 8.50e-11 @ Herbst et al. [1983] # R53 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + C2H4     -> C2H4+ + C             : 1.70e-10 @ Herbst et al. [1983] # R54 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + C2H4     -> C3H+  + H   +  H2     : 8.50e-11 @ Herbst et al. [1983] # R55 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + C2H4     -> C3H2+ + H2            : 3.40e-10 @ Herbst et al. [1983] # R56 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + C2H4     -> C3H3+ + H             : 1.00e-9 @ Herbst et al. [1983] # R57 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + C2H6     -> C2H2+ + CH4           : 1.70e-10 @ Herbst et al. [1983] # R58 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + C2H6     -> C2H3+ + CH3           : 5.10e-10 @ Herbst et al. [1983] # R59 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + C2H6     -> C2H5+ + CH            : 1.70e-10 @ Herbst et al. [1983] # R60 in Kim and Fox [1994] ")
reaction_rate_list.append(" C+    + C2H6     -> C3H3+ + H   +  H2     : 8.50e-10 @ Herbst et al. [1983] # R61 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH+   + H        -> C+    + H2            : 7.50e-10 @ Federer et al. [1984] # R62 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH+   + H2       -> CH2+  + H             : 1.20e-9 @ Federer et al. [1984] # R63 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH+   + CH4      -> C2H4+ + H             : 6.50e-11 @ Smith and Adams [1977] # R64 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH+   + CH4      -> C2H3+ + H2            : 1.10e-9 @ Smith and Adams [1977] # R65 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH+   + CH4      -> C2H2+ + H   +  H2     : 1.40e-10 @ Smith and Adams [1977] # R66 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH+   + C2H2     -> C3H2+ + H             : 2.40e-9 @ Anicich et al. [1986] # R67 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH2+  + H2       -> CH3+  + H             : 1.60e-9 @ Smith and Adams [1977] # R68 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH2+  + CH4      -> C2H5+ + H             : 3.60e-10 @ Smith and Adams [1977] # R69 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH2+  + CH4      -> C2H4+ + H2            : 8.40e-10 @ Smith and Adams [1977] # R70 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH2+  + C2H2     -> C3H3+ + H             : 2.50e-9 @ Anicich et al. [1986] # R71 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH3+  + CH4      -> C2H5+ + H2            : 1.20e-9 @ Adams and Smith [1978] # R72 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH3+  + C2H2     -> C3H3+ + H2            : 1.20e-9 @ Anicich et al. [1986] # R73 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH3+  + C2H4     -> C2H3+ + CH4           : 3.50e-10 @ Kim et al. [1977] # R74 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH3+  + C2H4     -> C3H3+ + H2  +  H2     : 4.60e-11 @ Kim et al. [1977] # R75 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH3+  + C2H4     -> C3H5+ + H2            : 5.20e-10 @ Kim et al. [1977] # R76 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH3+  + C2H6     -> C2H5+ + CH4           : 1.50e-9 @ Kim et al. [1977] # R77 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH3+  + C2H6     -> C3H5+ + H2  +  H2     : 1.60e-10 @ Kim et al. [1977] # R78 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH3+  + C2H6     -> C3H7+ + H2            : 1.00e-10 @ Kim et al. [1977] # R79 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH4+  + H2       -> CH5+  + H             : 3.0e-11 @ Federer et al. [1985] # R80 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH4+  + CH4      -> CH5+  + CH3           : 1.50e-9 @ Smith and Adams [1977] # R81 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH4+  + C2H2     -> C2H2+ + CH4           : 1.13e-9 @ Kim et al. [1977] # R82 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH4+  + C2H2     -> C2H3+ + CH3           : 1.23e-9 @ Kim et al. [1977] # R83 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH4+  + C2H2     -> C3H3+ + H   +  H2     : 1.51e-10 @ Kim et al. [1977] # R84 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH4+  + C2H4     -> C2H4+ + CH4           : 1.38e-9 @ Kim et al. [1977] # R85 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH4+  + C2H4     -> C2H5+ + CH3           : 4.23e-10 @ Kim et al. [1977] # R86 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH4+  + C2H4     -> C3H5+ + H   +  H2     : 5.52e-11 @ Kim et al. [1977] # R87 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH4+  + C2H6     -> C2H4+ + CH4 +  H2     : 1.91e-9 @ Kim et al. [1977] # R88 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH5+  + H        -> CH4+  + H2            : 1.50e-10 @ Federer et al. [1985] # R89 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH5+  + C2H2     -> C2H3+ + CH4           : 1.56e-9 @ Mackay et al. [1977] # R90 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH5+  + C2H4     -> C2H5+ + CH4           : 1.50e-9 @ Fiaux et al. [1974] # R91 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH5+  + C2H6     -> C2H5+ + CH4 +  H2     : 2.25e-10 @ Mackay et al. [1981] # R92 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH5+  + C2H6     -> C2H7+ + CH4           : 1.28e-9 @ Mackay et al. [1981] # R93 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2+   + H2       -> C2H+  + H             : 1.20e-9 @ Bohme and Wlodek [1990] # R94 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2+   + CH4      -> C2H+  + CH3           : 2.38e-10 @ Smith and Adams [1977] # R95 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2+   + CH4      -> C2H2+ + CH2           : 1.82e-10 @ Smith and Adams [1977] # R96 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2+   + CH4      -> C3H+  + H   +  H2     : 1.96e-10 @ Smith and Adams [1977] # R97 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2+   + CH4      -> C3H2+ + H2            : 5.74e-10 @ Smith and Adams [1977] # R98 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2+   + CH4      -> C3H3+ + H             : 2.10e-10 @ Smith and Adams [1977] # R99 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2+   + C2H2     -> C4H+  + H             : 1.20e-9 @ Knight et al. [1987] # R100 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2+   + C2H4     -> products              : 1.90e-9 @ Herod and Harrison [1970] # R101 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H+  + H2       -> C2H2+ + H             : 1.70e-9 @ Smith and Adams [1977] # R102 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H+  + CH4      -> C2H2+ + CH3           : 3.74e-10 @ Smith and Adams [1977] # R103 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H+  + CH4      -> C3H3+ + H2            : 3.74e-10 @ Smith and Adams [1977] # R104 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H+  + CH4      -> C3H4+ + H             : 1.32e-10 @ Smith and Adams [1977] # R105 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H+  + CH4      -> C3H5+                 : 2.20e-10 @ Smith and Adams [1977] # R106 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H+  + C2H2     -> C4H2+ + H             : 1.20e-9 @ Knight et al. [1987] # R107 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H+  + C2H4     -> products              : 1.71e-9 @ Herod and Harrison [1970] # R108 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + H2       -> C2H3+ + H             : 1.80e-12 @ Hansel et al. [1989] # R109 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + CH4      -> C3H4+ + H2            : 1.60e-10 @ Schiff et al. [1980] # R110 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + CH4      -> C3H5+ + H             : 6.40e-10 @ Schiff et al. [1980] # R111 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H2     -> C4H2+ + H2            : 5.16e-10 @ Knight et al. [1987] # R112 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H2     -> C4H3+ + H             : 5.88e-10 @ Knight et al. [1987] # R113 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H4     -> C2H4+ + C2H2          : 8.96e-10 @ Jarrold et al. [1983] # R114 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H4     -> C3H3+ + CH3           : 2.80e-10 @ Jarrold et al. [1983] # R115 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H4     -> C4H5+ + H             : 2.24e-10 @ Jarrold et al. [1983] # R116 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H6     -> C2H4+ + C2H4          : 2.63e-10 @ Kim et al. [1977] # R117 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H6     -> C2H5+ + C2H3          : 1.31e-10 @ Kim et al. [1977] # R118 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H6     -> C3H3+ + CH3 +  H2     : 8.76e-11 @ Kim et al. [1977] # R119 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H6     -> C3H4+ + CH4           : 1.46e-11 @ Kim et al. [1977] # R120 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H6     -> C3H5+ + CH3           : 7.88e-10 @ Kim et al. [1977] # R121 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H6     -> C4H5+ + H   +  H2     : 7.30e-11 @ Kim et al. [1977] # R122 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + C2H6     -> C4H7+ + H             : 1.31e-10 @ Kim et al. [1977] # R123 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H3+ + H        -> C2H2+ + H2            : 1.00e-10 @ Hansel et al. [1989] # R124 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H3+ + CH4      -> C3H5+ + H2            : 2.00e-10 @ Schiff et al. [1980] # R125 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H3+ + C2H2     -> C4H3+ + H2            : 2.16e-10 @ Anicich et al. [1986] # R126 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H3+ + C2H4     -> C2H5+ + C2H2          : 9.30e-10 @ Kim et al. [1977] # R127 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H3+ + C2H6     -> C2H5+ + C2H4          : 2.91e-10 @ Kim et al. [1977] # R128 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H3+ + C2H6     -> C3H5+ + CH4           : 2.48e-10 @ Kim et al. [1977] # R129 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H3+ + C2H6     -> C4H7+ + H2            : 8.06e-11 @ Kim et al. [1977] # R130 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H4+ + H        -> C2H3+ + H2            : 3.0e-10 @ Hansel et al. [1989] # R131 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H4+ + C2H2     -> C3H3+ + CH3           : 6.73e-10 @ Jarrold et al. [1983] # R132 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H4+ + C2H2     -> C4H5+ + H             : 2.37e-10 @ Jarrold et al. [1983] # R133 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H4+ + C2H4     -> C3H5+ + CH3           : 7.19e-10 @ Kim et al. [1977] # R134 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H4+ + C2H4     -> C4H7+ + H             : 7.11e-11 @ Kim et al. [1977] # R135 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H4+ + C2H6     -> C3H6+ + CH4           : 3.71e-13 @ Kim et al. [1977] # R136 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H4+ + C2H6     -> C3H7+ + CH3           : 4.93e-12 @ Kim et al. [1977] # R137 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H5+ + H        -> C2H4+ + H2            : 1.0e-11 @ Hansel et al. [1989] # R138 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H5+ + CH4      -> C3H7+ + H2            : 4.00e-15 @ Hiraoka and Kebarle [1975] # R139 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H5+ + C2H2     -> C3H3+ + CH4           : 6.84e-11 @ Kim et al. [1977] # R140 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H5+ + C2H2     -> C4H5+ + H2            : 1.22e-10 @ Kim et al. [1977] # R141 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H5+ + C2H4     -> C3H5+ + CH4           : 3.90e-10 @ Kim et al. [1977] # R142 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H5+ + C2H6     -> C4H9+ + H2            : 4.00e-11 @ Kim et al. [1977] # R143 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H6+ + H        -> C2H5+ + H2            : 1.0e-10 @ Hansel et al. [1989] # R144 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H6+ + C2H2     -> C2H5+ + C2H3          : 2.22e-10 @ Kim et al. [1977] # R145 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H6+ + C2H2     -> C3H5+ + CH3           : 8.19e-10 @ Kim et al. [1977] # R146 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H6+ + C2H2     -> C4H7+ + H             : 1.29e-10 @ Kim et al. [1977] # R147 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H6+ + C2H4     -> C2H4+ + C2H6          : 1.15e-9 @ Kim et al. [1977] # R148 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H6+ + C2H6     -> C3H8+ + CH4           : 7.98e-12 @ Kim et al. [1977] # R149 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H6+ + C2H6     -> C3H9+ + CH3           : 1.10e-11 @ Kim et al. [1977] # R150 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H7+ + C2H2     -> C2H3+ + C2H6          : 1.0e-9  @ estimate # R151 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H7+ + C2H4     -> C2H5+ + C2H6          : 1.0e-9  @ estimate # R152 in Kim and Fox [1994] ")

# Kim and Fox 1994 : 3 body reactions
reaction_rate_list.append(" H+    + H2 + H2 -> H3+   + H2 : 3.2e-29 @ Miller et al. [1968] # R153 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH3+  + H2 + H2 -> CH5+  + H2 : 3.3e-28 @ estimate # R154 in Kim and Fox [1994] ")
reaction_rate_list.append(" CH3+  + H2 + He -> CH5+  + He : 1.1e-28 @ Adams and Smith [1981] # R155 in Kim and Fox [1994] ")
reaction_rate_list.append(" C2H2+ + H2 + He -> C2H4+ + He : 1.2e-27 @ Adams and Smith [1981] # R156 in Kim and Fox [1994] ")
# Perry et al., 1999
reaction_rate_list.append(" H     + H  + H2 -> H2    + H2 : 2.7e-31 * (Tn/1)^(-0.6) @ Baulch et al. [1992] # Perry et al. [1999] ")

# Nakamura et al in prep
reaction_rate_list.append("  Na      +  hv      ->  Na+     +  e-             : Photoionization ")
reaction_rate_list.append("  Fe      +  hv      ->  Fe+     +  e-             : Photoionization ")
reaction_rate_list.append("  Mg      +  hv      ->  Mg+     +  e-             : Photoionization ")
reaction_rate_list.append("  Si      +  hv      ->  Si+     +  e-             : Photoionization ")

# Kim et al 2001
reaction_rate_list.append(" Meteoroid           ->  Na                        : Meteoroid ablation ")
reaction_rate_list.append(" Meteoroid           ->  Fe                        : Meteoroid ablation ")
reaction_rate_list.append(" Meteoroid           ->  Mg                        : Meteoroid ablation ")
reaction_rate_list.append(" Meteoroid           ->  Si                        : Meteoroid ablation ")
reaction_rate_list.append(" Meteoroid           ->  Na+     + e-              : Meteoroid ablation ")
reaction_rate_list.append(" Meteoroid           ->  Fe+     + e-              : Meteoroid ablation ")
reaction_rate_list.append(" Meteoroid           ->  Mg+     + e-              : Meteoroid ablation ")
reaction_rate_list.append(" Meteoroid           ->  Si+     + e-              : Meteoroid ablation ")

reaction_rate_list.append(" Na                  ->  Condensation              : 1.0e-5 ")
reaction_rate_list.append(" Fe                  ->  Condensation              : 1.0e-5 ")
reaction_rate_list.append(" Mg                  ->  Condensation              : 1.0e-5 ")
reaction_rate_list.append(" Si                  ->  Condensation              : 1.0e-5 ")

reaction_rate_list.append(" Na      +  H       +  H2      ->  NaH     +  H2  :  1.2e-30 * (Tn/1.0)^(-0.6) @ estimate")
reaction_rate_list.append(" Na      +  CH3     +  H2      ->  NaCH3   +  H2  :  4.0e-25 * (Tn/1.0)^(-1.8) @ estimate")
reaction_rate_list.append(" Na+     +  H2      +  H2      ->  NaH2+   +  H2  :  1.6e-30 * (Tn/80)^(-0.5) @ Smith et al. [1983]")
reaction_rate_list.append(" Na+     +  CH4     +  H2      ->  NaCH4+  +  H2  :  1.5e-28 * (Tn/80)^(-1.5) @ Smith et al. [1983]")
reaction_rate_list.append(" Na+     +  C2H2    +  H2      ->  NaC2H2+ +  H2  :  1.3e-28 * (Tn/80)^(-1.0) @ Smith et al. [1983]")
reaction_rate_list.append(" Na+     +  C2H4    +  H2      ->  NaC2H4+ +  H2  :  1.3e-28 * (Tn/80)^(-1.5) @ Smith et al. [1983]")
reaction_rate_list.append(" Na      +  H+      ->  Na+     +  H              :  1.2e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Na      +  H3+     ->  Na+     +  H       +  H2  :  1.1e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Na      +  CH3+    ->  Na+     +  CH3            :  3.4e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Na      +  C2H2+   ->  Na+     +  C2H2           :  2.7e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Na      +  C3H+    ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C3H2+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C3H3+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C3H4+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C3H5+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C3H6+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C3H7+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C3H8+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C3H9+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C4H+    ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C4H2+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C4H3+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C4H5+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C4H7+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  C4H9+   ->  Na+     +  products       :  3.0e-9 @ estimate")
reaction_rate_list.append(" Na      +  Si+     ->  Na+     +  Si             :  2.7e-9 @ Miller et al. [1997]")
#reaction_rate_list.append(" Na      +  H3O+    ->  Na+     +  H2O     +  H   :  3.1e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" NaH     +  H       ->  Na      +  H2             :  1.0e-10 * (Tn/300)^(-0.5) @ estimate")
reaction_rate_list.append(" Na+     +  e-      ->  Na                        :  2.7e-12 * (Te/300)^(-0.6) @ Miller et al. [1997]")
reaction_rate_list.append(" NaH2+   +  e-      ->  Na      +  H2             :  1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" NaH2+   +  e-      ->  NaH     +  H              :  1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" NaCH4+  +  e-      ->  Na      +  products       :  1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" NaC2H2+ +  e-      ->  Na      +  products       :  1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" NaC2H4+ +  e-      ->  Na      +  products       :  1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" NaCH4+  +  e-      ->  products                  :  1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" NaC2H2+ +  e-      ->  products                  :  1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" NaC2H4+ +  e-      ->  products                  :  1.5e-7  * (Te/300)^(-0.5) @ estimate")

reaction_rate_list.append(" Fe      +  H       +  H2      ->  FeH     +  H2  : 1.2e-30 * (Tn/1.0)^(-1.6) @ estimate")
reaction_rate_list.append(" Fe+     +  H       +  H2      ->  FeH+    +  H2  : 3.0e-30 * (Tn/1.0)^(-1.0) @ estimate")
reaction_rate_list.append(" Fe+     +  H2      +  H2      ->  FeH2+   +  H2  : 2.0e-30 @ Petrie et al. [1997]")
reaction_rate_list.append(" Fe+     +  CH4     +  H2      ->  FeCH4+  +  H2  : 1.1e-29 * (Tn/294)^(-1.5) @ Petrie et al. [1997]")
reaction_rate_list.append(" Fe+     +  C2H2    +  H2      ->  FeC2H2+ +  H2  : 1.4e-26 * (Tn/294)^(-1.0) @ Petrie et al. [1997]")
reaction_rate_list.append(" Fe+     +  C2H4    +  H2      ->  FeC2H4+ +  H2  : 1.4e-26 * (Tn/294)^(-1.5) @ Petrie et al. [1997]")
reaction_rate_list.append(" Fe+     +  C2H6    ->  Fe      +  C2H4+   +  H2  : 1.1e-12 * (Tn/300)^(-0.6) @ Schultz et al. [1988]")
reaction_rate_list.append(" Fe+     +  Na      ->  Fe      +  Na+            : 1.0e-11 @ Miller et al. [1997]")
reaction_rate_list.append(" FeH+    +  H       ->  Fe+     +  H2             : 1.8e-9 @ estimate")
reaction_rate_list.append(" FeH2+   +  H       ->  FeH+    +  H2             : 1.8e-9 @ estimate")
reaction_rate_list.append(" Fe      +  H+      ->  Fe+     +  H              : 7.4e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Fe      +  H3+     ->  Fe+     +  H       +  H2  : 4.9e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Fe      +  C+      ->  Fe+     +  C              : 2.6e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Fe      +  CH3+    ->  Fe+     +  CH3            : 2.4e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Fe      +  C2H2+   ->  Fe+     +  C2H2           : 2.0e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Fe      +  C3H+    ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C3H2+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C3H3+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C3H4+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C3H5+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C3H6+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C3H7+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C3H8+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C3H9+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C4H+    ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C4H2+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C4H3+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C4H5+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C4H7+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  C4H9+   ->  Fe+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Fe      +  Si+     ->  Fe+     +  Si             : 1.9e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" FeH     +  H       ->  Fe      +  H2             : 1.0e-10 * (Tn/300)^( 0.5) @ estimate")
reaction_rate_list.append(" Fe+     +  e-      ->  Fe                        : 3.7e-12 * (Te/300)^(-0.6) @ Miller et al. [1997]")
reaction_rate_list.append(" FeH+    +  e-      ->  Fe      +  H              : 3.0e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" FeH2+   +  e-      ->  FeH     +  H              : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" FeH2+   +  e-      ->  Fe      +  H2             : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" FeCH4+  +  e-      ->  Fe      +  products       : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" FeC2H2+ +  e-      ->  Fe      +  products       : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" FeC2H4+ +  e-      ->  Fe      +  products       : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" FeCH4+  +  e-      ->  products                  : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" FeC2H2+ +  e-      ->  products                  : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" FeC2H4+ +  e-      ->  products                  : 1.5e-7  * (Te/300)^(-0.5) @ estimate")

reaction_rate_list.append(" Mg      +  H       +  H2      ->  MgH     +  H2  : 1.2e-30 * (Tn/1.0)^(-0.6) @ estimate")
reaction_rate_list.append(" Mg+     +  H       +  H2      ->  MgH+    +  H2  : 3.0e-30 * (Tn/1.0)^(-1.0) @ Smith et al. [1983]")
reaction_rate_list.append(" Mg+     +  H2      +  H2      ->  MgH2+   +  H2  : 1.6e-30 * (Tn/80.0)^(-0.5) @ Smith et al. [1983]")
reaction_rate_list.append(" Mg+     +  CH4     +  H2      ->  MgCH4+  +  H2  : 1.5e-28 * (Tn/80.0)^(-1.5) @ Smith et al. [1983]")
reaction_rate_list.append(" Mg+     +  C2H2    +  H2      ->  MgC2H2+ +  H2  : 1.3e-28 * (Tn/80.0)^(-1.0) @ Smith et al. [1983]")
reaction_rate_list.append(" Mg+     +  C2H4    +  H2      ->  MgC2H4+ +  H2  : 1.3e-28 * (Tn/80.0)^(-1.5) @ Smith et al. [1983]")
reaction_rate_list.append(" Mg+     +  Na      ->  Mg      +  Na+            : 1.0e-11 @ Miller et al. [1997]")
reaction_rate_list.append(" MgH+    +  H       ->  Mg+     +  H2             : 1.8e-9 @ estimate")
reaction_rate_list.append(" MgH2+   +  H       ->  MgH+    +  H2             : 1.8e-9 @ estimate")
reaction_rate_list.append(" Mg      +  H+      ->  Mg+     +  H              : 1.1e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Mg      +  H3+     ->  Mg+     +  H       +  H2  : 1.0e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Mg      +  C+      ->  Mg+     +  C              : 1.1e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Mg      +  CH3+    ->  Mg+     +  CH3            : 3.5e-10 @ Miller et al. [1997]")
reaction_rate_list.append(" Mg      +  CH5+    ->  Mg+     +  CH4     +  H   : 1.4e-9 @ Po and Porter [1977]")
reaction_rate_list.append(" Mg      +  C2H2+   ->  Mg+     +  C2H2           : 3.0e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" Mg      +  C3H+    ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C3H2+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C3H3+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C3H4+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C3H5+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C3H6+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C3H7+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C3H8+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C3H9+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C4H+    ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C4H2+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C4H3+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C4H5+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C4H7+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  C4H9+   ->  Mg+     +  products       : 3.0e-9 @ estimate")
reaction_rate_list.append(" Mg      +  Si+     ->  Mg+     +  Si             : 2.9e-9 @ Miller et al. [1997]")
reaction_rate_list.append(" MgH     +  H       ->  Mg      +  H2             : 1.0e-10 * (Tn/300)^( 0.5) @ estimate")
reaction_rate_list.append(" Mg+     +  e-      ->  Mg                        : 2.8e-12 * (Te/300)^(-0.8) @ Miller et al. [1997]")
reaction_rate_list.append(" MgH+    +  e-      ->  Mg      +  H              : 3.0e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" MgH2+   +  e-      ->  Mg      +  H2             : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" MgH2+   +  e-      ->  MgH     +  H              : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" MgCH4+  +  e-      ->  Mg      +  products       : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" MgC2H2+ +  e-      ->  Mg      +  products       : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" MgC2H4+ +  e-      ->  Mg      +  products       : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" MgCH4+  +  e-      ->  products                  : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" MgC2H2+ +  e-      ->  products                  : 1.5e-7  * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" MgC2H4+ +  e-      ->  products                  : 1.5e-7  * (Te/300)^(-0.5) @ estimate")

reaction_rate_list.append(" Si      +  H       +  H2      ->  SiH     +  H2  : 1.2e-30 * (Tn/1.0)^(-0.6) @ estimate")
reaction_rate_list.append(" Si+     +  H       +  H2      ->  SiH+    +  H2  : 3.0e-30 * (Tn/1.0)^(-1.0) @ estimate")
reaction_rate_list.append(" Si+     +  H2      +  H2      ->  SiH2+   +  H2  : 1.6e-30 * (Tn/80)^(-0.5) @ Smith et al. [1983]")
reaction_rate_list.append(" Si+     +  CH4     +  H2      ->  SiCnHm+ +  H2  : 1.4e-28 * (Tn/80)^(-1.5) @ Smith et al. [1983]") # <- modified:actually Si+ + CH4 + H2 -> SiCH4+ + H2
reaction_rate_list.append(" Si+     +  C2H2    +  H2      ->  SiCnHm+ +  H2  : 1.2e-28 * (Tn/80)^(-1.0) @ Smith et al. [1983]") # <- modified:actually Si+ + C2H2 + H2 -> SiC2H2+ + H2
reaction_rate_list.append(" Si+     +  C2H4    +  H2      ->  SiCnHm+ +  H2  : 1.2e-28 * (Tn/80)^(-1.5) @ Smith et al. [1983]") # <- modified:actually Si+ + C2H4 + H2 -> SiC2H4+ + H2
reaction_rate_list.append(" Si+     +  H2      ->  SiH2+                      : 3.0e-18 @ Herbst et al. [1989]")
reaction_rate_list.append(" Si+     +  CH3     ->  SiCnHm+ +  H               : 1.0e-9  @ Herbst et al. [1989]") # <- modified:actually Si+ + CH3 -> SiCH2+ + H
reaction_rate_list.append(" Si+     +  CH3     ->  SiCnHm+ +  H2              : 1.0e-9  @ Herbst et al. [1989]") # <- modified:actually Si+ + CH3 -> SiCH+ + H2
reaction_rate_list.append(" Si+     +  CH4     ->  SiCnHm+ +  H               : 7.7e-11 @ Cheng et al. [1973]") # <- modified:actually Si+ + CH4 -> SiCH3+ + H
reaction_rate_list.append(" Si+     +  C2H2    ->  SiCnHm+ +  H               : 1.8e-10 @ Wlodek et al. [1991]") # <- modified:actually Si+ + C2H2 -> SiC2H+ + H
reaction_rate_list.append(" Si+     +  C2H4    ->  SiCnHm+ +  CH3     +  H    : 7.4e-11 @ Wlodek et al. [1991]") # <- modified:actually Si+ + C2H4 -> SiC+ + CH3 + H
reaction_rate_list.append(" Si+     +  C2H4    ->  SiCnHm+ +  H               : 1.0e-9  @ Herbst et al. [1989]") # <- modified:actually Si+ + C2H4 -> SiC2H3+ + H
reaction_rate_list.append(" Si+     +  C2H6    ->  SiCnHm+ +  CH4             : 1.2e-10 @ Wlodek et al. [1991]") # <- modified:actually Si+ + C2H6 -> SiCH2+ + CH4
reaction_rate_list.append(" Si+     +  C2H6    ->  SiCnHm+ +  CH3             : 6.4e-10 @ Wlodek et al. [1991]") # <- modified:actually Si+ + C2H6 -> SiCH3+ + CH3
reaction_rate_list.append(" Si+     +  C2H6    ->  SiCnHm+ +  CH3     +  H2   : 2.4e-11 @ Wlodek et al. [1991]") # <- modified:actually Si+ + C2H6 -> SiCH2+ + CH3 + H2
#reaction_rate_list.append(" Si+     +  H2O     ->  SiOH+   +  H               : 2.3e-10 @ Fahey et al. [1981]")
reaction_rate_list.append(" SiH+    +  H       ->  Si+     +  H2              : 1.9e-9  @ Herbst et al. [1989]")
reaction_rate_list.append(" SiH2+   +  H       ->  SiH+    +  H2              : 1.8e-9  @ estimate")
reaction_rate_list.append(" Si      +  H+      ->  Si+     +  H               : 9.9e-10 @ Miller et al. [1997]")
reaction_rate_list.append(" Si      +  H3+     ->  SiH+    +  H2              : 3.7e-9  @ Miller et al. [1997]")
reaction_rate_list.append(" Si      +  CH3+    ->  SiCnHm+ +  H               : 2.0e-9  @ Herbst et al. [1989]") # <- modified:actually Si + CH3+ -> SiCH2+ + H
reaction_rate_list.append(" Si      +  CH3+    ->  SiCnHm+ +  H2              : 2.0e-9  @ Herbst et al. [1989]") # <- modified:actually Si + CH3+ -> SiCH+ + H2
reaction_rate_list.append(" Si      +  CH5+    ->  SiCnHm+ +  H               : 2.0e-10 @ Herbst et al. [1989]") # <- modified:actually Si + CH5+ -> SiCH4+ + H
reaction_rate_list.append(" Si      +  CH5+    ->  SiCnHm+ +  H2              : 2.0e-10 @ Herbst et al. [1989]") # <- modified:actually Si + CH5+ -> SiCH3+ + H2
reaction_rate_list.append(" Si      +  C2H3+   ->  SiCnHm+ +  H2              : 2.0e-10 @ Herbst et al. [1989]") # <- modified:actually Si + C2H3+ -> SiC2H+ + H2
reaction_rate_list.append(" Si      +  C2H3+   ->  SiCnHm+ +  H               : 2.0e-10 @ Herbst et al. [1989]") # <- modified:actually Si + C2H3+ -> SiC2H2+ + H
reaction_rate_list.append(" Si      +  C3H+    ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C3H2+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C3H3+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C3H4+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C3H5+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C3H6+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C3H7+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C3H8+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C3H9+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C4H+    ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C4H2+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C4H3+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C4H5+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C4H7+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
reaction_rate_list.append(" Si      +  C4H9+   ->  SiCnHm+ +  products        : 2.0e-10 @ estimate")
#reaction_rate_list.append(" Si      +  H3O+    ->  SiH+    +  H3O             : 1.8e-9  @ Miller et al. [1997]")
reaction_rate_list.append(" Si      +  CH2     ->  SiCH    +  H               : 2.0e-11 * (Tn/300)^(-0.5) @ Herbst et al. [1989]")
reaction_rate_list.append(" Si      +  CH3     ->  SiCH2   +  H               : 4.0e-11 * (Tn/300)^( 0.5) @ Herbst et al. [1989]")
reaction_rate_list.append(" SiH     +  H       ->  Si      +  H2              : 1.0e-10 * (Tn/300)^( 0.5) @ estimate")
reaction_rate_list.append(" Si+     +  e-      ->  Si                         : 4.9e-12 * (Te/300)^(-0.6) @ Miller et al. [1997]")
reaction_rate_list.append(" SiH+    +  e-      ->  Si      +  H               : 2.0e-7 * (Te/300)^(-0.5) @ Miller et al. [1997]")
reaction_rate_list.append(" SiH2+   +  e-      ->  Si      +  H       +  H    : 2.0e-7 * (Te/300)^(-0.5) @ Miller et al. [1997]")
reaction_rate_list.append(" SiH2+   +  e-      ->  SiH     +  H               : 1.5e-7 * (Te/300)^(-0.5) @ Miller et al. [1997]")
reaction_rate_list.append(" SiH2+   +  e-      ->  Si      +  H2              : 1.5e-7 * (Te/300)^(-0.5) @ Miller et al. [1997]")
reaction_rate_list.append(" SiCnHm+ +  e-      ->  Si                         : 1.5e-7 * (Te/300)^(-0.5) @ estimate")
reaction_rate_list.append(" SiCnHm+ +  e-      ->  products                   : 1.5e-7 * (Te/300)^(-0.5) @ estimate")



#--------------------------------------------------------------------------------------------------------------------------
#
#                                                       Saturn
#
#--------------------------------------------------------------------------------------------------------------------------
#Planet_list.append(['Saturn',len(reaction_rate_list)])
# reaction list




# --------------------------------------------------------------------------------------------------------------------------
#
#                                                        Titan
#
# --------------------------------------------------------------------------------------------------------------------------
#Planet_list.append(['Titan',len(reaction_rate_list)])
# reaction list



# --------------------------------------------------------------------------------------------------------------------------
#
#                                               END of Reaction list
#
# --------------------------------------------------------------------------------------------------------------------------
Planet_list.append(['END',len(reaction_rate_list)])
# --------------------------------------------------------------------------------------------------------------------------



###########################################################################################################################
#
#                                                  Reference List: Label
#
###########################################################################################################################
doi_list.append(["Anicich [1993a]", "Anicich, V. G. (1993), Evaluated bimolecular ion-molecule gas-phase kinetics of positive-ions for use in modeling planetary-atmospheres, cometary comae, and interstellar clouds, J. Phys. Chem. Ref. Data, 22, 1469-1569, doi:10.1063/1.555940."])
doi_list.append(["Anicich [1993b]", "Anicich, V. G. (1993), A survey of bimolecular ion-molecule reactions for use in modeling the chemistry of planetary atmospheres, cometary comae, and interstellar clouds, J. Suppl. Ser., 84, 215-315, doi:10.1086/191752."])
doi_list.append(["Chaffin et al. [2017]", "Chaffin, M. S., J. Deighan, N. M. Schneider, and A. I. F. Stewart (2017), Elevated atmospheric escape of atomic hydrogen from Mars induced by high-altitude water, Nature Geoscience, 10, 174-178, doi:10.1038/ngeo2887."])
doi_list.append(["Fox and Sung [2001]", "Fox, J. L., and Sung, K. Y. (2001), Solar activity variations of the Venus thermosphere/ionosphere, J. Geophys. Res., 106(  A10), 21305-21335, doi:10.1029/2001JA000069."])
doi_list.append(["Kim and Fox [1994]", "Kim, Y. H., and J. L. Fox (1994), The chemistry of hydrocarbon ions in the Jovian ionosphere, Icarus, 112, 310-325, doi:10.1006/icar.1994.1186."])
doi_list.append(['Kim et al. [2001]', 'Kim, Y. H., W. D. Pesnell, J. M. Grebowsky, and J. L. Fox (2001), Meteoric Ions in the Ionosphere of Jupiter, Icarus, 150, 261-278, doi:10.1006/icar.2001.6590.'])
doi_list.append(["Molina-Cuberos et al. [2001]", "Molina-Cuberos, G. J., J. J. López-Moreno, R. Rodrigo, H. Lichtenegger, and K. Schwingenschuh (2001). A model of the Martian ionosphere below 70 km, Adv. Space Res., 27(11), 1801-1806, doi:10.1016/S0273-1177(01)00342-8."])
doi_list.append(["Mukundan et al. [2020]", "Mukundan, V., S. V. Thampi, A. Bhardwaj, and C. Krishnaprasad (2020), The dayside ionosphere of Mars: Comparing a one-dimensional photochemical model with MAVEN Deep Dip campaign observations, Icarus, 337, 113502, doi:10.1016/j.icarus.2019.113502."])
doi_list.append(["Perry et al. [1999]", "Perry, J. J., Y. H. Kim, J. L. Fox, and H. S. Porter (1999), Chemistry of the Jovian auroral ionosphere, J. Geophys. Res., 104(E7), 16541-16565, doi:10.1029/1999JE900022."])
doi_list.append(["Pinto et al. [1980]", "Pinto, J. P., G. R. Gladstone, and Y. L. Yung (1980). Photochemical production of formaldehyde in Earth's primitive atmosphere, Science, 210, 183-185, doi:10.1126/science.210.4466.183."])
doi_list.append(["Schunk and Nagy [1980]", "Schunk, R. W., and A. F. Nagy (1980), Ionospheres of the terrestrial planets, Rev. Geophys., 18, 813-852, doi:10.1029/RG018i004p00813."])
doi_list.append(["Schunk and Nagy [2000]", "Schunk, R. W., and A. F. Nagy (2000), Ionospheres: Physics, Plasma Physics, and Chemistry, Cambridge University Press, doi:10.1017/CBO9780511635342."])
doi_list.append(['Schunk and Nagy [2000]', 'Schunk, R. W., and A. F. Nagy (2000), Ionization and energy exchange processes, in Ionosphere Physics, Plasma Physics, and Chemistry, edited by J. T. Houghton, M. J. Rycroft, and A. J. Dessler, pp. 237-267 and pp. 521-531, Cambridge Univ. Press, Cambridge, U. K.'])
doi_list.append(["Shinagawa et al. [1987]", "Shinagawa, H., T. E. Cravens, and A. F. Nagy (1987), A one-dimensional time-dependent model of the magnetized ionosphere of Venus, J. Geophys. Res., 92, 7317-7330, doi:10.1029/JA092iA07p07317."])
doi_list.append(["Verronen et al. [2016]", "Verronen, P. T., M. E. Andersson, D. R. Marsh, T. Kovács, and J. M. C. Plane (2016), WACCM-D—Whole Atmosphere Community Climate Model with D-region ion chemistry, J. Adv. Model. Earth Syst., 8, 954-975, doi:10.1002/2015MS000592."])


###########################################################################################################################
#
#                                                Reference List: Experiments
#
###########################################################################################################################
doi_list.append(["Adams and Smith [1981]", "Adams, N. G., and D. Smith (1981), The rate coefficients for several ternary association reactions at thermal energies., Chem. Phys. Lett., 79, 563-567, doi:10.1016/0009-2614(81)85036-1."])
doi_list.append(["Adams and Smith [1987]", "Adams, N. G., and D. Smith (1987), Recent advances in the studies of reaction rates relevant to interstellar chemistry, Astrochemistry IAU Symposium N. 120 (M. Vardya and S. Tarafdar, Eds.) pp. 1-18, Reidel, Dordrecht"])
doi_list.append(["Adams and Smith [1988]", "Adams, N. G., and D. Smith (1988), Measuremets of the dissociative recombination coefficients for several polyatomic ion species at 300 K, Chem. Phys. Lett., 144, 11-14, doi:10.1016/0009-2614(88)87081-7."])
doi_list.append(["Adams et al. [1984]", "Adams, N. G., D. Smith, and E. Alge (1984), Measurements of dissociative recombination coefficients of H3+, HCO+, N2H+, and CH5+ at 95 and 300 K using the FALP apparatus, J. Chem. Phys., 81, 1778-1784, doi:10.1063/1.447849."])
doi_list.append(["Alge et al. [1983]", "Alge, E., N. G. Adams, and D. Smith (1983), Measurements of the dissociative recombination coefficients of O2+, NO+ and NH4+ in the temperature range 200-600 K, J. Phys. B, At. Mol. Phys., 16, 1433-1444, doi:10.1088/0022-3700/16/8/017."])
doi_list.append(["Anicich and Huntless [1986]", "Anicich, V. G., and W. T. Huntress (1986), A survey of bimolecular ion-molecule reactions for use in modeling the chemistry of planetary atmospheres, cometary comae, and interstellar clouds, Astrophys. J. Suppl. Ser., 62, 553-672, doi:10.1086/191151."])
doi_list.append(["Anicich et al. [1986]", "Anicich, V. G., W. T. Huntress, Jr., and M. J. McGowan (1986), Ion-molecule reactions of hydrocarbon ions in acetylene and hydrocyanic acid, J. Phys. Chem., 90, 2446-2450, doi:10.1021/j100402a038."])
doi_list.append(["Atkinson and Welge [1972]", "Atkinson, R., and K. H. Welge (1972), Temperature dependence of O(1S) deactivation by CO2, O2, N2, and Ar, J. Chem. Phys., 57, 3689-3693, doi:10.1063/1.1678829."])
doi_list.append(["Atkinson et al. [1997]", "Atkinson, R., D. L. Bauch, R. A. Cox, R. F. Hampson Jr., L. A. Kerr, M. J. Rossi, and J. Troe (1997), Evaluated kinetic and photochemical data for atmospheric chemistry, suppl. VI, J. Phys. Chem. Ref. Data, 26, 1329-1499, doi:10.1063/1.556010."])
doi_list.append(["Auerbach et al. [1977]", "Auerbach, D., R. Cacak, R. Candano, C. Keyser, T. Gaily, J. B. A. Mitchell, J. Wm. McGowan, and S. F. J. Wilk (1977), Merged electron-ion beam experiments. 1. Method and Measurements of (e-H2+) and (e-H3+) dissociative recombination cross sections, J. Phys., B10, 3797-3820, doi:10.1088/0022-3700/10/18/033."])
doi_list.append(["Bates [1989]", "Bates, D. R. (1989), Theoretical considerations regarding some inelastic atomic collision processes of interest in aeronomy : Deactivation and charge transfer, Planet. Space Sci., 17, 363-368, doi:10.1016/0032-0633(89)90033-0."])
doi_list.append(["Bates and Dalgarno [1962]", "Bates, D. R., and A. Dalgarno (1962), Electron recombination, in Atomic Molicular Processes, pp. 245-271, Academic Press, New York, doi:10.1016/B978-0-12-081450-3.50011-4."])
doi_list.append(["Baulch et al. [1992]", "Baulch, D. L., C. J. Cobos, R. A. Cox, C. Esser, P. Frank, T. Just, J. A. Kerr, M. J. Pilling, J. Troe, R. W. Walker, J. Warnatz (1992), Evaluated kinetic data for combustion modeling, J. Phys. Chem. Ref. Data, 21, 411-734, doi:10.1063/1.555908."])
doi_list.append(["Baulch et al. [1994]", "Baulch, D. L., C. J. Cobos, R. A. Cox, P. Frank, G. Hayman, T. Just, J. A. Kerr, T. Murrells, M. J. Pilling, T. Troe, R. W. Walker, J. Warnatz (1994), Evaluated kinetic data for combustion modeling supplement-1, J. Phys. Chem. Ref. Data, 23, 847-1033, doi:10.1063/1.555953."])
doi_list.append(["Berrington and Burke [1981]", "Berrington, K. A., and P. G. Burke (1981), Effective collision strengths for forbideen transitions in e- - N and e- - O scattering, Planet. Space Sci.,29, 377-381, doi:10.1016/0032-0633(81)90026-X."])
doi_list.append(["Bischof and Linder [1986]", "Bischof, G., and F. Linder (1986), Crossed beam study of He+ - O2 charge transfer reactions in the collision energy range 0.5-200 eV, Z. Phys. D, 1, 303-320, doi:10.1007/BF01436687."])
doi_list.append(["Bohme and Wlodek [1990]", "Bohme, D. K., and S. Wlodek (1990), Hydrogenation of carbon-cluster cations with molecular hydrogen: implications for the growth of carbon cluster molecules, Int. J. Mass Spectrom. Ion. Proc., 102, 133-149 , doi:10.1016/0168-1176(90)80056-9."])
doi_list.append(["Bohme et al. [1980]", "Bohme, D. K., and G. I. Mackay, and H. I. Schiff (1980), Determination of proton assinities from the kinetics of proton tansfer reactions. VII. The proton affinities of O2, H2, Kr, O, N2, Xe, CO2, CH4, N2O, and CO, J. Chem. Phys., 73, 4976-4986, doi:10.1063/1.439975."])
doi_list.append(["Bohme et al. [1982]", "Bohme, D. K., A. B. Rakshit, and H. I. Schiff (1982), Reactions of 12C+ with hydrocarbons at 296 K: Carbon-carbon bond formation, Chem. Phys. Lett. 93, 592-597, doi:10.1016/0009-2614(82)83736-6."])
doi_list.append(["Canosa et al. [1992]", "Canosa, A., J. C. Gomet, B. R. Rowe, J. B. A. Mitchell, and J. L. Queffelec (1992), Furthermeasurements of the H3+(v = 0, 1, 2) dissociative recombination rate coefficient, J. Chem. Phys., 97, 1028-1037, doi:10.1063/1.463282."])
doi_list.append(["Capetanakis et al. [1993]", "Capetanakis, F. P., F. Sonderman, S. Hoser, and F. Stuhl (1993), Temperature dependence of the quenching of O(1S) by simple inorganic molecules, J. Geophys. Res., 98, 7883-3887 , doi:10.1063/1.464596."])
doi_list.append(["Chastaing et al. [2000]", "Chastaing, D. , S. D. Le Picard, and I. R. Sims (2000), Direct kinetic measurements on reactions of atomic carbon, C(3P), with O2 and NO at temperatures down to 15 K, J. Chem. Phys., 112, 8466-8469, doi:10.1063/1.481448."])
doi_list.append(["Clow and Futrell [1972]", "Clow, R. P., and J. H. Futrell (1972), Ion-molecule reactions in isotopic hydrogen by ion cyclotron resonance, Int. J. Mass Spetctrom. Ion Phys., 10, 405-418, doi:10.1016/0020-7381(72)80003-2."])
doi_list.append(["Constantinides et al. [1979]", "Constantinides, E. R., J. H. Black, and A. Dalgarno (1979), The photochemistry of N+ ions, Geophys. Res. Lett., 61, 569-572, doi:10.1029/GL006i007p00569."])
doi_list.append(["Dalgarno [1979]", "Dalgarno, A. (1979), A. E. reaction rate data, Rep. AFGL-TR 79-0067, Air Force Geophys. Lab., Bedford, Mass"])
doi_list.append(["Dheandhanoo et al. [1984]", "Dheandhanoo, S., R. Johnson, and M. A. Biondi (1984), Measured ion-molecule reaction rates for modelling Titan's atmosphere, Planet. Space Sci., 32, 1301-1305, doi:10.1016/0032-0633(84)90073-4."])
doi_list.append(["Dotan et al. [1997]", "Dotan, I., P. M. Hierl, R. A. Morris, and A. A. Viggiano (1997), Rate constants for the reactions of N+ and N2+ with O2 as a function of temperature (300 - 1800 K), Int. J. Mass Spectrom. Ion Proc., 167/168, 223-230, doi:10.1016/S0168-1176(97)00077-3."])
doi_list.append(["Dotan et al. [1999]", "Dotan, I., A. J. Midey, and A. A. Viggiano (1999), Rate constants for the reactions of Ar+ with CO2 and SO2 as a function of themperature (300-1500 K), J. Am. Soc. Mass Spectrom., 10, 815-820, doi:10.1016/S1044-0305(99)00056-2."])
doi_list.append(["Dotan et al. [2000]", "Dotan, I., A. J. Midey, and A. A. Viggiano (2000), Kinetics of the reactions of N2+ with CO2 and SO2 from 300-1400 K, J. Chem. Phys., 113, 1732-1737, doi:10.1063/1.481975."])
doi_list.append(["Dotan and Lindinger [1982]", "Dotan, I. and W. Lindinger (1982), Energy dependences of thereactions of Ar+ with H2, N2, CO, O2, CO2, N2O, and COS, J. Chem. Phys., 76, 4972-4977, doi:10.1063/1.442843."])
doi_list.append(["Fahey et al. [1981a]", "Fahey, D. W., I. Dotan, F. C. Fehesenfeld, D. L. Albritton, and L. A. Viehland (1981a), Energy dependense of the rate constant of the reaction N+ + NO at collision energies 0.04 to 2.5 eV, J. Chem. Phys., 74, 3320-3323, doi:10.1063/1.441484."])
doi_list.append(["Fahey et al. [1981b]", "Fahey, D. W., F. C. Fehsenfeld, and E. E. Ferguson (1981b), Rate constant for the reaction C+ + CO2 at collision energies 0.04 to 2.5 eV, Geophys. Res. Lett., 8, 1115-1117, doi:10.1029/GL008i010p01115."])
doi_list.append(["Federer et al. [1984]", "Federer, W., H. Villinger, F. Howorka, W. Lindinger, P. Tosi, D. Bassi, and E. Ferguson (1984), Reaction of O+, CO+, and CH+ ions with atomic hydrogen, Phys. Rev. Lett., 52, 2084-2086 , doi:10.1103/PhysRevLett.52.2084."])
doi_list.append(["Federer et al. [1985]", "Federer, W., H. Villinger, P. Tosi, D. Bassi, E. E. Ferguson, and W. Lindinger (1985), Laboratory studies of ion reactions with atomic hydrogen, Molecular Astrophysics (G. H. F. Diercksen et al. Ed.), pp. 649-655, Reidel, Boston"])
doi_list.append(["Fehsenfeld and Ferguson [1972]", "Fersenfeld, F. C., and E. E. Ferguson (1972), Thermal energy reaction rates constants for H+ and CO+ with O and NO, J. Chem. Phys., 56, 3066-3070, doi:10.1063/1.1677642."])
doi_list.append(["Fehsenfeld et al. [1970]", "Fersenfeld, F. C., D. B. Dunkin, and E. E. Ferguson (1970), Rate constants for the reaction of CO2+ with O, O2, and NO; N2+ with O and NO; and O2+ with NO, Planet. Space Sci., 18, 1267-1269, doi:10.1016/0032-0633(70)90216-3."])
doi_list.append(["Fell et al. [1990]", "Fell, C., J. Stienfeld, and S. Miller (1990), Quenching of N(2D) by O(3P), J. Chem. Phys., 92, 4768-4777, doi:10.1063/1.457694."])
doi_list.append(["Ferguson [1973]", "Ferguson, E. E. (1973), Rate constants of thermal energy binary ion-molecule reactions of aeronomic interest, At. Data Nucl. Data Tables, 12, 159-178, doi:10.1016/0092-640X(73)90017-X."])
doi_list.append(["Ferguson et al. [1992]", "Ferguson, E. E., J. M. Van Doren, A. A. Viggino, and R. A. Morris (1992), Internal and tarnslational energy effects on the charge-transfer reaction of CO2+ with O2, Int. J. Mass Spectrom. Ion Proc., 117, 261-282, doi:10.1016/0168-1176(92)80098-L."])
doi_list.append(["Fiaux et al. [1974]", "Fiaux, A., D. L. Smith, and J. H. Futrell (1974), Reaction of CH5+ with C2H2, C2H4, C3H6 and c-C3H6, Int. J. Mass Spectrom. Ion Phys., 15, 9-21, doi:10.1016/0020-7381(74)80082-3."])
doi_list.append(["Fox [1982a]", "Fox, J. L. (1982), The chemistry of metastable species in the Venusian ionosphere, Icarus, 51, 248-260, doi:10.1016/0019-1035(82)90081-1."])
doi_list.append(["Fox [2015]", "Fox, J. L. (2015), The chemistry of protonated species in the Martian ionosphere, Icarus, 252, 366-392, doi:10.1016/j.icarus.2015.01.010."])
doi_list.append(["Froese-Fischer and Saha [1983]", "Froese-Fischer, C., and H. P. Saha (1983), Multiconfiguration Hartree-Fock results with Breit-Pauli corrections for forbidden transitions in the 2P4 configuration, Phys. Rev. A, 28, 3169-3178, doi:10.1103/PhysRevA.28.3169."])
doi_list.append(["Frost et al. [1998]", "Frost, M. J., S. Kato, V. M. Bierbaum, and S. R. Leone (1998), Reactions of N2+(v) with CO and NO at thermal energy, Chem. Phys. 231, 145-153, doi:10.1016/S0301-0104(97)00348-0."])
doi_list.append(["Gerlich [1992]", "Gerlich, D. (1992), Inhomogeneous RF fields: A versatile tool for the study of processes with slow ions, in State-Selected and State-to-State Ion-Molecule Reaction Dynamics, Part 1, Experiment, edited by C.-Y. Ng and M. Baer, Adv. Chem. Phys. Ser., vol. LXXXII, John Wiley, New York, doi:10.1002/9780470141397.ch1."])
doi_list.append(["Glosik et al. [1978]", "Glosik, J., A. B. Rakshit, N. D. Twiddy, N. G. Adams, and D. Smith (1978), Measurement of the rates of reaction of the ground and metastable excited states of O2+, NO+ and O+ with atmospheric gases at thermal energy, J. Phys. B, At. Mol. Phys., 11, 3365-3379 , doi:10.1088/0022-3700/11/19/013."])
doi_list.append(["Goldan et al. [1966]", "Goldan, P. D., A. L. Schmeltekopf, F. C. Fehsenfeld, H. I. Schiff, and E. E. Ferguson (1966), Thermal energy ion-molecule reaction rates, 2, Some reactions of ionospheric interest, J. Chem. Phys., 44, 4095-4013, doi:10.1063/1.1726588."])
doi_list.append(["Gougousi et al. [1997]", "Gougousi, T., M. F. Golde, and R. Johnsen (1997), Electron-ion recombination rate coefficient measurements in a flowing afterglow plasma, Chem. Phys. Lett., 265, 399-403, doi:10.1016/S0009-2614(96)01488-1."])
doi_list.append(["Hansel et al. [1989]", "Hansel, A., R. Richter, W. Lindinger, and E. E. Ferguson (1989), Reactions of C2 and C3 hydrocarbon ions with H, D, H2 and D2 at near-thermal energies, Int. J. Mass Spectrom. Ion Proc, 94, 251-260, doi:10.1016/0168-1176(89)80046-1."])
doi_list.append(["Herbst et al. [1983]", "Herbst, E., N. G. Adams, and D. Smith (1983), Laboratory measurements of ion-molecular reactions pertaining to interstellar hydrocarbon synthesis, Astrophys. J., 269, 329-333, doi:10.1086/161046."])
doi_list.append(["Herod and Harrison [1970]", "Herod, A. A., and A. G. Harrison (1970), Bimolecular reactions of ions trapped in an electron space charge, Int. J. Mass Spectrom. Ion Phys. 4, 415-431, doi:10.1016/0020-7381(70)85054-9."])
doi_list.append(["Herron [1999]", "Herron, J. T. (1999), Evaluated chemical kinetic data for reactions of N(2D), N(2P), and N2(A3**) in the gas phase, J. Phys. Chem. Ref. Data, 28, 1453-1483, doi:10.1063/1.556043."])
doi_list.append(["Hierl et al. [1997]", "Hierl, P. M., I. Dotan, J. V. Seeley, J. M. Van Doren, R. A. Morris, and A. A. Viggiano (1997), Rate constants for the reactions of O+ with N2 and O2 as a function of temeprature (300 - 1800 K), J. Chem. Phys., 106, 3540-3544, doi:10.1063/1.473450."])
doi_list.append(["Hiraoka and Kebarle [1975]", "Hiraoka, K., and P. Kebarle (1975), Temperature dependence of bimolecular ion molecule reactions. The reaction C2H5+ + CH4 = C3H7+ + H2, J. Chem. Phys., 63, 394-397, doi:10.1063/1.431116."])
doi_list.append(["Huntress et al. [1974]", "Huntress, W. T., J. K. Kim, and L. P. Thread (1974), On the reaction of protons with methane, ammonia, water and oxygen, Chem. Phys. Lett., 29, 189-190, doi:10.1016/0009-2614(74)85009-8."])
doi_list.append(["Jamieson et al. [1992]", "Jamieson, M. J., J. M. Finch, R. S. Friedman, and A. Dalgarno (1992), Collisional excitation of metastable oxygen O(1D) atoms through the 3\u03A3+u channel of O2, Planet. Space Sci., 40, 1719-1721, doi:10.1016/0032-0633(92)90128-B."])
doi_list.append(["Jarrold et al. [1983]", "Jarrold, M. F., L. M. Bass, P. R. Kemper, P. A. M. Van Koppen, and M. T. Bowers (1983), Unimolecular and bimolecular reactions in the C4H6+ system: Experiment and theory, J. Chem. Phys., 78, 3756-3766, doi:10.1063/1.445151."])
doi_list.append(["Johnsen and Biondi [1980]", "Johnsen, R., and M. A. Biondi (1980), Laboratory measurements of the O+(2D) + N2 and O+(2D) + O2 reaction rate coefficients and their ionospheric implications, Gesphys. Res. Lett., 7, 401-403, doi:10.1029/GL007i005p00401."])
doi_list.append(["Jones et al. [1981]", "Jones, J. D. C., K. Birkinshaw, and N. D. Twiddy (1981), Rate coefficients and product ion distributions for the reactions of OH+ and H2O+ with  N2, O2, NO, N2O, Xe, CO, CO2, H2S and H2 at 300 K, Chem. Phys. Lett., 77, 484-488, doi:10.1016/0009-2614(81)85191-3."])
doi_list.append(["Jones et al. [1986]", "Jones, M. E., S. E. Barlow, G. B. Ellison, and E. E. Ferguson (1986), Reactions of C+, He+, and Ne+ with vibrationally excited H2 and D2, Chem. Phys. Lett. 130, 218-223, doi:10.1016/0009-2614(86)80458-4."])
doi_list.append(["Karpas et al. [1979]", "Karpas, Z., V. G. Anicich, and W. T. Huntress (1979), An ion cylotron resonance study of reactions of ions with hydrogen atoms, J. Chem. Phys., 70, 2877-2881, doi:10.1063/1.437823."])
doi_list.append(["Kella et al. [1996]", "Kella, D., P. J. Johnson, H. B. Pedersen, L. VejbyChristensen, and L. H. Andersen (1996), Branching ratios for dissociative recombination of 15N14N, Phys. Rev. Lett., 77, 2432-2435, doi:10.1103/PhysRevLett.77.2432."])
doi_list.append(["Kella et al. [1997]", "Kella, D., L. Vejby-Christenson, P. J. Johnson, H. B. Pedersen, and L. H. Andersen (1997), The source of green light emission determined from a heavy-ion storage ring experiment, Science, 276, 1530-1533, doi:10.1126/science.276.5318.1530."])
doi_list.append(["Kernahan and Pang [1975]", "Kernahan, J. A., and P. H. L. Pang (1975), Experimental determination of absolute A coefficients for 'forbidden' atomi oxygen lines, Can. J. Phys., 53, 455-458, 1975, doi:10.1139/p75-058."])
doi_list.append(["Kim and Huntress [1975a]", "Kim, J. K., and W. T. Huntress (1975), Product distribution and rate constants for the reactions of thermal energy He+ ions with some neutral hydrides and hydrocarbons, Int. J. Mass. Spectrom. Ion Phys., 16, 451-454, doi:10.1016/0020-7381(75)85034-0."])
doi_list.append(["Kim and Huntress [1975b]", "Kim, J. K., and W. T. Huntress (1975), Ion cyclotron resonance studies on the reaction of H2+ and D2+ ions with various simple molecules and hydrocarbons, J. CHem. Phys., 62, 2820-2825, doi:10.1063/1.430817."])
doi_list.append(["Kim et al. [1977]", "Kim, J. K., V. G. Anicich, and W. T. Huntress (1977), Product distribution and rate constants for the reactions of CH3+, CH4+, C2H2+, C2H3+, C2H4+, C2H5+, and C2H6+ ions with CH4, C2H2, C2H4, and C2H6, J. Phys. Chem., 81, 1798-1805, doi:10.1021/j100534a002."])
doi_list.append(["Knight et al. [1987]", "Knight, J. S., C. G. Freeman, M. J. McEwan, V. G. Anicich, and W. T. Huntress, Jr. (1987), A flowtube study of ion-molecule reactions of acetylene, J. Phys. Chem., 91, 3898-3902, doi:10.1021/j100298a033."])
doi_list.append(["Lee et al. [1978]", "Lee, J. H., J. V. Michael, W. A. Payne, and L. Stief (1978), Absolute rate of the reaction of N(4S) with NO from 196-400 K with DF-RF and FP-RF techniques, J. Chem. Phys., 69, 3069-3076, doi:10.1063/1.436998."])
doi_list.append(["Le-Teuff et al. [1999]", "Le-Teuff, Y. H., T. J. Millar, and A. J. Markwick (1999), The UMIST database for astrochemistry, Astron. Astrophys. Suppl., 146, 157-168, doi:10.1051/aas:2000265."])
doi_list.append(["Li et al. [1997a]", "Li, X., Y.-L. Huang, G. D. Flesch, and C. Y. Ng (1997), Absolute state-selected total cross sections for the ion-molecule reactions O+(4S, 2D, 2P) + H2(D2), J. Chem. Phys., 106, 564-571, doi:10.1063/1.473395."])
doi_list.append(["Li et al. [1997b]", "Li, X., Y.-L. Huang, G. D. Flesch, and C. Y. Ng (1997), A state-selected study of the ion-molecule reactions O+(4S, 2D, 2P) + N2, J. Chem. Phys., 106,1373-1381, doi:10.1063/1.474087."])
doi_list.append(["Mackay et al. [1977]", "Mackay, G. I., K. Tanaka, and D. K. Bohme (1977), Rate constants at 297 K for proton-transfer reactions with acetylen: An assessment of the average quadrupole orientation theory, Int. J. Mass Spectrom. Ion Phys., 24, 125-136, doi:10.1016/0020-7381(77)80020-X."])
doi_list.append(["Mackay et al. [1981]", "Mackay, G. I., H. I. Schiff, and D. K. Bohme (1981), A room-temperature study of the kinetic and energetics for the protonation of ethane, Can. J. Chem., 59, 1771-1778, doi:10.1139/v81-265."])
doi_list.append(["Matta et al. [2013]", "Matta, M., Withers, P., and M. Mendillo (2013), The composition of mars' topside ionosphere: Effects of hydrogen, J. Geophys. Res., 118,2693-5681, doi:10.1002/jgra.50104."])
doi_list.append(["Mauclaire et al. [1978]", "Mauclaire, G., R. Derai, and R. Marx (1978), ICR determination of kinetically excited ions produced in water and methane by charge transfer from thermal rare gas ions, Int. J. Mass Spectrom. Ion Phys., 26, 289-301, doi:10.1016/0020-7381(78)80031-X."])
doi_list.append(["McElroy and McConnell [1971]", "McElroy, M. B., and M. C. McConnell, Atomic carbon in the atmospheres of Mars and Venus, J. Geophys. Res., 76, 6674-6690, doi:10.1029/JA076i028p06674."])
doi_list.append(["McFarland et al. [1974]", "McFarland, M., D. L. Albritton, F. C. Fehesenfeld, E. E. Ferguson, and A. L. Schmeltekopf (1974), Energy dependence and branching ratios of the N2+ + O reaction, J. Geophys. Res., 79, 2925-2926, doi:10.1029/JA079i019p02925."])
doi_list.append(["McLaughlin and Bell [1998]", "McLaughhn, B. M., and K. L. Bell (1998), Electron-impact excitation of the fine-structure levels (1s2 2s2 2p3 4PO3/2, 2DO5/2,3/2, 2PO3/2,1/2) of singly ionized atomic oxygen, J. Phys. B, At. Mol. Phys., 31, 4317-4329, doi:10.1088/0953-4075/31/19/017."])
doi_list.append(["Mehr and Biondi [1969]", "Mehr, F. J., and M. A. Biondi, Electron temperature dependence of recombination of O2+ and N2+ ions with electrons, Phys. Rev., 181, 264-270, doi:10.1103/PhysRev.181.264."])
doi_list.append(["Midey and Viggiano [1998]", "Midey, A. J., and A. A. Viggiano (1998), Rate constants for the reaction of Ar+ with O2 and CO as function of temperature from 300 K to 1400 K: Derivation of rotational and vibrational energy effects, J. Chem. Phys., 109, 5257-5263, doi:10.1063/1.477142."])
doi_list.append(["Midey and Viggiano [1999]", "Midey, A. J., and A. A. Viggiano(1999), Rate constants for the reaction of O2+ with NO from 300 to 1400 K, J. Chem. Phys., 110, 10,746-10,748, doi:10.1063/1.479017."])
doi_list.append(["Miller et al. [1968]", "Miller, T. M., J. T. Moseley, D. W. Martine, and E. W. McDaniel (1968), Reactions of H+ in H2 and D+ in D2: Mobilities of hydrogen and alkali ions in H2 and D2 gases, Phys. Rev., 173, 115-123, doi:10.1103/PhysRev.173.115."])
doi_list.append(["Miller et al. [1984b]", "Miller, T. M., R. E. Wetterskog, and J. F. Paulson (1984), Temperature dependence of the ion-molecule reactions N+ + CO, C+ + NO, and C+, CO+, CO2+ + O2 from 90-450K, J. Chem. Phys., 80, 4922-4925, doi:10.1063/1.446514."])
doi_list.append(["Mitchell and McGowan [1978]", "Mitchell, J. B. A., and J. Wm. McGowan (1978), The dissociative recombination of CH+ X1\u03A3  + (v = 0), Astrophys. J. Lett., 222, L77-L79, doi:10.1086/182696."])
doi_list.append(["Mitchell et al. [1983]", "Mitchell, J. B. A., J. L. Forand, C. T. Ng, D. P. Levac, R. E. Mitchell, P. M. Mul, W. Claeye, A. Sen, and J. Wm. McGowan (1983), Measurements of the branching ratio for the dissociative recombination of H3+ + e, Phys. Rev. Lett., 51, 885-888, doi:10.1103/PhysRevLett.51.885."])
doi_list.append(["Mul and McGowan [1980]", "Mul, P. M., and J. Wm. McGowan (1980), Dissociative recombination of C2+, C2H+, C2H2+, and C2H3+, Astrophys. J., 237, 749-751, doi:10.1086/157921."])
doi_list.append(["Mul et al. [1981]", "Mul, P. M., J. B. A. Mitchell, V. S. D'Angelo, P. Defrance, J. Wm. McGowan, and H. R. Froelich (1981), Merged electron-ion bean experiments. IV. Dissociative recombination for the methane group CH+, CH2+, CH3+, CH4+, CH5+, J. Phys., B 14, 1353-1361, doi:10.1088/0022-3700/14/8/020."])
doi_list.append(["O'Keefe et al. [1986]", "O'Keefe, A., G. Mauclaire, D. Parent, and M. T. Bowers (1986), Product energy disposal in the reaction of N+(3P) with O2(X3 \u03A3), J. Chem. Phys., 84, 215-219, doi:10.1063/1.450173."])
doi_list.append(["Piper [1993]", "Piper, L. G. (1993), The reactions of N(2P) with O2 and O, J. Chem. Phys., 98, 8560-8564, doi:10.1063/1.464515."])
doi_list.append(["Prasad and Huntress [1980]", "Prasad, S. S., and W. T. Huntress (1980), A modelfor gas phase chemistry in interstellar clouds: 1. The basic model, library of chemical reactions, and chemistry among C, N and O compounds, Astrophys. J. Suppl. Ser., 43, 1-35, doi:10.1086/190665."])
doi_list.append(["Rakshit [1982]", "Rakshit, A. B. (1982), A drift-chamber mass-spectrometric study of the interaction of triatomic (1+) ions with netral molecules at 300 K, Int. J. Mass Spectrom. Ion Phys., 41, 185-197, doi:10.1016/0020-7381(82)85034-1."])
doi_list.append(["Rawlins et al. [1989]", "Rawlins, W. T., M. E. Fraser, and S. M. Miller (1989), Robvidrational excitation of nitric oxide in the reaction of O2 with atomic nitrogen, J. Phys, Chem., 93, 1097-1107, doi:10.1021/j100340a016."])
doi_list.append(["Rosen et al. [1998]", "Rosen, S., R. Peverall, M. Larsson, A. LePadellec, J. Semaniak, A. Larson, C. Stromholm, W. J. van der Zande, H. Danared, and G. H. Dunn (1998), Absolute cross sections and final-state distributions for dissociative recombination and excitation of CO+(v-0) using an ion storage ring, Phys. Rev. A, 57, 4462-4471, doi:10.1103/PhysRevA.57.4462."])
doi_list.append(["Rusch et al. [1977]", "Rusch, D. W., D. G. Torr, P. B. Hays, and J. C. G. Walker (1977), The OII (7319-7330 A) dayglow, J. Geophys. Res., 82, 719-722, doi:10.1029/JA082i004p00719."])
doi_list.append(["Schauer et al. [1989]", "Schauer, M. M., S. R. Jefferts, S. E. Barlow, and G. H. Dunn (1989), Reactions of H2 with He+ at temperature below40 K, J. Chem. Phys., 91, 4593-4596, doi:10.1063/1.456748."])
doi_list.append(["Schiff et al. [1980]", "Schiff, H. I., G. I. Mackay, G. D. Vlachos, and D. K. Bohme (1980), Interstallar Molecules (B. H. Andrew Ed.), pp. 307-310, D. Reidel, Boston"])
doi_list.append(["Schofield [1978]", "Schofield, K. (1978), Rate constants for the gaseous interactions of O(21D2) and O(21S0): A critical evaluation, J. Photochem., 9, 55-68, doi:10.1016/0047-2670(78)87006-3."])
doi_list.append(["Scott et al. [1997]", "Scott, G. B. I., D. A. Fairley, C. G. Freeman, M. J. McEwan, P. Spanel, and D. Smith (1997), Gas phase reactions of some positice ions with atomic and molecular hydrogen at 300 K, J. Chem. Phys., 106, 3982-3987, doi:10.1063/1.473116."])
doi_list.append(["Scott et al. [1998]", "Scott, G. B. I., D. A. Fairley, C. G. Freeman, M. J. McEwan, and V. J. Anicich (1998), Gas-phasereactions of some positie ions with atomic and molecular nitrogen, J. Chem. Phys., 109, 9010-9014, doi:10.1063/1.477571."])
doi_list.append(["Scott et al. [1999]", "Scott, G. B. I., D. A. Fairley, D. B. Milligan, C. G. Freeman, and M. J. McEwan (1999), Gas phase reactions of some positive ions with atomic and molecular oxygen and nitric oxide at 300 K, J. Phys. Chem. A, 103, 7470-7473, doi:10.1021/jp9913719."])
doi_list.append(["Seaton and Osterbrock [1957]", "Seaton, M. J., and D. E. Osterbrock (1957), Relative intensities in gaseous nebulae, Astrophys. J., 125, 66-83, 1957, doi:10.1086/146282."])
doi_list.append(["Shihira et al. [1994]", "Shihira, Y., T. Suzuki, S. Unayama, H. Umemoto, and S. Tsunashima (1994), Reactions of N(22D) and N(22P) with O2, J. Chem. Soc. Faraday  Trans., 90, 549-552, doi:10.1039/ft9949000549."])
doi_list.append(["Smith and Adams [1977]", "Smith, D., and N. G. Adams (1977), Reaction of simple hydrocarbon ions with molecules at thermal energies, Int. J. Mass Spectrom. Ion Phys., 23, 123-135, doi:10.1016/0020-7381(77)80094-6."])
doi_list.append(["Smith and Futrell [1976]", "Smith, R. D., and J. H. Futrell (1976), Reactions of rare-gas hydride ions with ethylene, Int. J. Mass Spectrom. Ion Phys., 20, 71-76, doi:10.1016/0020-7381(76)80033-2."])
doi_list.append(["Smith et al. [1976]", "Smith, R. D., D. L. Smith, and J. H. Futrell (1976), Ion-molecule reactions in H2/rare-gas sysems by ion cyclotron resonance. 1. Reactions with He and Ne, Int. J. Mass Spectrom. Ion Phys., 19, 369-394, doi:10.1016/0020-7381(76)80020-4."])
doi_list.append(["Streit et al. [1976]", "Streit, G. E., C. J. Howard, and A. L. Schmeltekopf (1976), Temperature dependence of O(1D) rate constants for reactions with O2, N2, CO, O3 and H2O, J. Chem. Phys., 65, 4761-4764, doi:10.1063/1.432930."])
doi_list.append(["Torr and Torr [1980]", "Torr, D. G., and M. R. Torr (1980), Determination of the thermal rate coefficient, products and branching ratios for the reaction of O+(2D) with N2, J. Geophys. Res., 85, 783-787, doi:10.1029/JA085iA02p00783."])
doi_list.append(["Tsang and Hampson [1986]", "Tsang, W., and R. F. Hampson, Chemical kinetic data base for combustion chemistry, part 1, Methane and related compounds, J. Phys. Chem. Ref. Data, 15, 1087-1279, doi:10.1063/1.555759."])
doi_list.append(["Uiterwaal et al. [1995]", "Uiterwaal, C. J. G. J., J. van Eck, and A. Niehaus (1995), State-selected ion-molecule reactions: Charge transfer and atomic rearrangement processes in thermal energy collisions of H2+(X;v) + N2 and of N2+(X, A;v) + H2, J. Chem. Phys., 102, 744-753, doi:10.1063/1.469187."])
doi_list.append(["Vejby-Christensen et al. [1998]", "Vejby-Christensen, L., D. Kella, H. B. Pedersen, and L. H. Andersen (1998), Dissociative recombination of NO+, Phys. Rev. A, 57, 3627-3634, doi:10.1103/PhysRevA.57.3627."])
doi_list.append(["Viggiano et al. [1990]", "Viggiano, A. A., R. A. Morris, and J. F. Paulson (1990), Rate constant and branching fractions for the reaction of O+(2D, 2P) with CO2, J. Chem. Phys., 93, 1483-1484, doi:10.1063/1.459163."])
doi_list.append(["Wiese et al. [1966]", "Wiese, W. L., M. W. Smith, and B. M. Glennon (1966), Atomic transition probabilities, Ref. Data Ser. Natl. Bur. Stand., U.S., 4"])
doi_list.append(["Young and Dunn [1975]", "Young, R. A., and O. J. Dunn (1975), Excitation and quenching of N(2P), J. Chem. Phys., 63, 1150-1153, doi:10.1063/1.431441."])
doi_list.append(["Yousif and Mitchell [1989]", "Yousif, F. B., and J. B. A. Mitchell (1989), Recombination and excitation of HeH+, Phys. Rev., A 40, 4318-4321, doi:10.1103/PhysRevA.40.4318."])
doi_list.append(["Zipf [1980]", "Zipf, E. C. (1980), The dissociative recombination of vibrationally excited N2+ ions, Geophys. Res. Lett., 7, 645-648, doi:10.1029/GL007i009p00645."])




###########################################################################################################################
# reaction -> unicode for display
# visionally make everything good by unicode
def upper_upper_unicode(char):
    # ^**
    char = re.sub('\u20700' , '\u2070\u2070', char)
    char = re.sub('\u20701' , '\u2070\u00B9', char)
    char = re.sub('\u20702' , '\u2070\u00B2', char)
    char = re.sub('\u20703' , '\u2070\u00B3', char)
    char = re.sub('\u20704' , '\u2070\u2074', char)
    char = re.sub('\u20705' , '\u2070\u2075', char)
    char = re.sub('\u20706' , '\u2070\u2076', char)
    char = re.sub('\u20707' , '\u2070\u2077', char)
    char = re.sub('\u20708' , '\u2070\u2078', char)
    char = re.sub('\u20709' , '\u2070\u2079', char)

    char = re.sub('\u00B90' , '\u00B9\u2070', char)
    char = re.sub('\u00B91' , '\u00B9\u00B9', char)
    char = re.sub('\u00B92' , '\u00B9\u00B2', char)
    char = re.sub('\u00B93' , '\u00B9\u00B3', char)
    char = re.sub('\u00B94' , '\u00B9\u2074', char)
    char = re.sub('\u00B95' , '\u00B9\u2075', char)
    char = re.sub('\u00B96' , '\u00B9\u2076', char)
    char = re.sub('\u00B97' , '\u00B9\u2077', char)
    char = re.sub('\u00B98' , '\u00B9\u2078', char)
    char = re.sub('\u00B99' , '\u00B9\u2079', char)

    char = re.sub('\u00B20' , '\u00B2\u2070', char)
    char = re.sub('\u00B21' , '\u00B2\u00B9', char)
    char = re.sub('\u00B22' , '\u00B2\u00B2', char)
    char = re.sub('\u00B23' , '\u00B2\u00B3', char)
    char = re.sub('\u00B24' , '\u00B2\u2074', char)
    char = re.sub('\u00B25' , '\u00B2\u2075', char)
    char = re.sub('\u00B26' , '\u00B2\u2076', char)
    char = re.sub('\u00B27' , '\u00B2\u2077', char)
    char = re.sub('\u00B28' , '\u00B2\u2078', char)
    char = re.sub('\u00B29' , '\u00B2\u2079', char)

    char = re.sub('\u00B30' , '\u00B3\u2070', char)
    char = re.sub('\u00B31' , '\u00B3\u00B9', char)
    char = re.sub('\u00B32' , '\u00B3\u00B2', char)
    char = re.sub('\u00B33' , '\u00B3\u00B3', char)
    char = re.sub('\u00B34' , '\u00B3\u2074', char)
    char = re.sub('\u00B35' , '\u00B3\u2075', char)
    char = re.sub('\u00B36' , '\u00B3\u2076', char)
    char = re.sub('\u00B37' , '\u00B3\u2077', char)
    char = re.sub('\u00B38' , '\u00B3\u2078', char)
    char = re.sub('\u00B39' , '\u00B3\u2079', char)

    char = re.sub('\u20740' , '\u2074\u2070', char)
    char = re.sub('\u20741' , '\u2074\u00B9', char)
    char = re.sub('\u20742' , '\u2074\u00B2', char)
    char = re.sub('\u20743' , '\u2074\u00B3', char)
    char = re.sub('\u20744' , '\u2074\u2074', char)
    char = re.sub('\u20745' , '\u2074\u2075', char)
    char = re.sub('\u20746' , '\u2074\u2076', char)
    char = re.sub('\u20747' , '\u2074\u2077', char)
    char = re.sub('\u20748' , '\u2074\u2078', char)
    char = re.sub('\u20749' , '\u2074\u2079', char)

    char = re.sub('\u20750' , '\u2075\u2070', char)
    char = re.sub('\u20751' , '\u2075\u00B9', char)
    char = re.sub('\u20752' , '\u2075\u00B2', char)
    char = re.sub('\u20753' , '\u2075\u00B3', char)
    char = re.sub('\u20754' , '\u2075\u2074', char)
    char = re.sub('\u20755' , '\u2075\u2075', char)
    char = re.sub('\u20756' , '\u2075\u2076', char)
    char = re.sub('\u20757' , '\u2075\u2077', char)
    char = re.sub('\u20758' , '\u2075\u2078', char)
    char = re.sub('\u20759' , '\u2075\u2079', char)

    char = re.sub('\u20760' , '\u2076\u2070', char)
    char = re.sub('\u20761' , '\u2076\u00B9', char)
    char = re.sub('\u20762' , '\u2076\u00B2', char)
    char = re.sub('\u20763' , '\u2076\u00B3', char)
    char = re.sub('\u20764' , '\u2076\u2074', char)
    char = re.sub('\u20765' , '\u2076\u2075', char)
    char = re.sub('\u20766' , '\u2076\u2076', char)
    char = re.sub('\u20767' , '\u2076\u2077', char)
    char = re.sub('\u20768' , '\u2076\u2078', char)
    char = re.sub('\u20769' , '\u2076\u2079', char)

    char = re.sub('\u20770' , '\u2077\u2070', char)
    char = re.sub('\u20771' , '\u2077\u00B9', char)
    char = re.sub('\u20772' , '\u2077\u00B2', char)
    char = re.sub('\u20773' , '\u2077\u00B3', char)
    char = re.sub('\u20774' , '\u2077\u2074', char)
    char = re.sub('\u20775' , '\u2077\u2075', char)
    char = re.sub('\u20776' , '\u2077\u2076', char)
    char = re.sub('\u20777' , '\u2077\u2077', char)
    char = re.sub('\u20778' , '\u2077\u2078', char)
    char = re.sub('\u20779' , '\u2077\u2079', char)

    char = re.sub('\u20780' , '\u2078\u2070', char)
    char = re.sub('\u20781' , '\u2078\u00B9', char)
    char = re.sub('\u20782' , '\u2078\u00B2', char)
    char = re.sub('\u20783' , '\u2078\u00B3', char)
    char = re.sub('\u20784' , '\u2078\u2074', char)
    char = re.sub('\u20785' , '\u2078\u2075', char)
    char = re.sub('\u20786' , '\u2078\u2076', char)
    char = re.sub('\u20787' , '\u2078\u2077', char)
    char = re.sub('\u20788' , '\u2078\u2078', char)
    char = re.sub('\u20789' , '\u2078\u2079', char)

    char = re.sub('\u20790' , '\u2079\u2070', char)
    char = re.sub('\u20791' , '\u2079\u00B9', char)
    char = re.sub('\u20792' , '\u2079\u00B2', char)
    char = re.sub('\u20793' , '\u2079\u00B3', char)
    char = re.sub('\u20794' , '\u2079\u2074', char)
    char = re.sub('\u20795' , '\u2079\u2075', char)
    char = re.sub('\u20796' , '\u2079\u2076', char)
    char = re.sub('\u20797' , '\u2079\u2077', char)
    char = re.sub('\u20798' , '\u2079\u2078', char)
    char = re.sub('\u20799' , '\u2079\u2079', char)

    return char

def lower_lower_unicode(char):
    # ^**
    char = re.sub('\u20800' , '\u2080\u2080', char)
    char = re.sub('\u20801' , '\u2080\u2081', char)
    char = re.sub('\u20802' , '\u2080\u2082', char)
    char = re.sub('\u20803' , '\u2080\u2083', char)
    char = re.sub('\u20804' , '\u2080\u2084', char)
    char = re.sub('\u20805' , '\u2080\u2085', char)
    char = re.sub('\u20806' , '\u2080\u2086', char)
    char = re.sub('\u20807' , '\u2080\u2087', char)
    char = re.sub('\u20808' , '\u2080\u2088', char)
    char = re.sub('\u20809' , '\u2080\u2089', char)

    char = re.sub('\u20810' , '\u2081\u2080', char)
    char = re.sub('\u20811' , '\u2081\u2081', char)
    char = re.sub('\u20812' , '\u2081\u2082', char)
    char = re.sub('\u20813' , '\u2081\u2083', char)
    char = re.sub('\u20814' , '\u2081\u2084', char)
    char = re.sub('\u20815' , '\u2081\u2085', char)
    char = re.sub('\u20816' , '\u2081\u2086', char)
    char = re.sub('\u20817' , '\u2081\u2087', char)
    char = re.sub('\u20818' , '\u2081\u2088', char)
    char = re.sub('\u20819' , '\u2081\u2089', char)

    char = re.sub('\u20820' , '\u2082\u2080', char)
    char = re.sub('\u20821' , '\u2082\u2081', char)
    char = re.sub('\u20822' , '\u2082\u2082', char)
    char = re.sub('\u20823' , '\u2082\u2083', char)
    char = re.sub('\u20824' , '\u2082\u2084', char)
    char = re.sub('\u20825' , '\u2082\u2085', char)
    char = re.sub('\u20826' , '\u2082\u2086', char)
    char = re.sub('\u20827' , '\u2082\u2087', char)
    char = re.sub('\u20828' , '\u2082\u2088', char)
    char = re.sub('\u20829' , '\u2082\u2089', char)

    char = re.sub('\u20830' , '\u2083\u2080', char)
    char = re.sub('\u20831' , '\u2083\u2081', char)
    char = re.sub('\u20832' , '\u2083\u2082', char)
    char = re.sub('\u20833' , '\u2083\u2083', char)
    char = re.sub('\u20834' , '\u2083\u2084', char)
    char = re.sub('\u20835' , '\u2083\u2085', char)
    char = re.sub('\u20836' , '\u2083\u2086', char)
    char = re.sub('\u20837' , '\u2083\u2087', char)
    char = re.sub('\u20838' , '\u2083\u2088', char)
    char = re.sub('\u20839' , '\u2083\u2089', char)

    char = re.sub('\u20840' , '\u2084\u2080', char)
    char = re.sub('\u20841' , '\u2084\u2081', char)
    char = re.sub('\u20842' , '\u2084\u2082', char)
    char = re.sub('\u20843' , '\u2084\u2083', char)
    char = re.sub('\u20844' , '\u2084\u2084', char)
    char = re.sub('\u20845' , '\u2084\u2085', char)
    char = re.sub('\u20846' , '\u2084\u2086', char)
    char = re.sub('\u20847' , '\u2084\u2087', char)
    char = re.sub('\u20848' , '\u2084\u2088', char)
    char = re.sub('\u20849' , '\u2084\u2089', char)

    char = re.sub('\u20850' , '\u2085\u2080', char)
    char = re.sub('\u20851' , '\u2085\u2081', char)
    char = re.sub('\u20852' , '\u2085\u2082', char)
    char = re.sub('\u20853' , '\u2085\u2083', char)
    char = re.sub('\u20854' , '\u2085\u2084', char)
    char = re.sub('\u20855' , '\u2085\u2085', char)
    char = re.sub('\u20856' , '\u2085\u2086', char)
    char = re.sub('\u20857' , '\u2085\u2087', char)
    char = re.sub('\u20858' , '\u2085\u2088', char)
    char = re.sub('\u20859' , '\u2085\u2089', char)

    char = re.sub('\u20860' , '\u2086\u2080', char)
    char = re.sub('\u20861' , '\u2086\u2081', char)
    char = re.sub('\u20862' , '\u2086\u2082', char)
    char = re.sub('\u20863' , '\u2086\u2083', char)
    char = re.sub('\u20864' , '\u2086\u2084', char)
    char = re.sub('\u20865' , '\u2086\u2085', char)
    char = re.sub('\u20866' , '\u2086\u2086', char)
    char = re.sub('\u20867' , '\u2086\u2087', char)
    char = re.sub('\u20868' , '\u2086\u2088', char)
    char = re.sub('\u20869' , '\u2086\u2089', char)

    char = re.sub('\u20870' , '\u2087\u2080', char)
    char = re.sub('\u20871' , '\u2087\u2081', char)
    char = re.sub('\u20872' , '\u2087\u2082', char)
    char = re.sub('\u20873' , '\u2087\u2083', char)
    char = re.sub('\u20874' , '\u2087\u2084', char)
    char = re.sub('\u20875' , '\u2087\u2085', char)
    char = re.sub('\u20876' , '\u2087\u2086', char)
    char = re.sub('\u20877' , '\u2087\u2087', char)
    char = re.sub('\u20878' , '\u2087\u2088', char)
    char = re.sub('\u20879' , '\u2087\u2089', char)

    char = re.sub('\u20880' , '\u2088\u2080', char)
    char = re.sub('\u20881' , '\u2088\u2081', char)
    char = re.sub('\u20882' , '\u2088\u2082', char)
    char = re.sub('\u20883' , '\u2088\u2083', char)
    char = re.sub('\u20884' , '\u2088\u2084', char)
    char = re.sub('\u20885' , '\u2088\u2085', char)
    char = re.sub('\u20886' , '\u2088\u2086', char)
    char = re.sub('\u20887' , '\u2088\u2087', char)
    char = re.sub('\u20888' , '\u2088\u2088', char)
    char = re.sub('\u20889' , '\u2088\u2089', char)

    char = re.sub('\u20890' , '\u2089\u2080', char)
    char = re.sub('\u20891' , '\u2089\u2081', char)
    char = re.sub('\u20892' , '\u2089\u2082', char)
    char = re.sub('\u20893' , '\u2089\u2083', char)
    char = re.sub('\u20894' , '\u2089\u2084', char)
    char = re.sub('\u20895' , '\u2089\u2085', char)
    char = re.sub('\u20896' , '\u2089\u2086', char)
    char = re.sub('\u20897' , '\u2089\u2087', char)
    char = re.sub('\u20898' , '\u2089\u2088', char)
    char = re.sub('\u20899' , '\u2089\u2089', char)

    return char

def reaction_unicode(char):
    nlabel0 = 0 # ex) if " ... nH2O" -> nlabel0 = 1
    char = re.sub('->' , ' \u2192 ', char)
    if re.compile('\s\d').search(char):
        num0 = re.findall('\s+(\d{1,3})', char)
        nlabel0 = len(num0)
        char = re.sub('\s+\d{1,3}', ' ?na?', char,1)
        char = re.sub('\s+\d{1,3}', ' ?nb?', char,1)
        char = re.sub('\s+\d{1,3}', ' ?nc?', char,1)
        char = re.sub('\s+\d{1,3}', ' ?nd?', char,1)
        char = re.sub('\s+\d{1,3}', ' ?ne?', char,1)
        char = re.sub('\s+\d{1,3}', ' ?nf?', char,1)
        char = re.sub('\s+\d{1,3}', ' ?ng?', char,1)
        char = re.sub('\s+\d{1,3}', ' ?nh?', char,1)
        char = re.sub('\s+\d{1,3}', ' ?ni?', char,1)
        char = re.sub('\s+\d{1,3}', ' ?nj?', char,1)
        #print(num0)

    nlabel1 = 0 # ex) if "nH2O" -> nlabel1 = 1
    if re.compile('^\d').search(char):
        num1 = re.findall('^\d{1,3}', char)
        nlabel1 = 1
        char = re.sub('^\d{1,3}', '?ma?', char,1)
        print(num1, char)

    # For isotope expression ex) ^13C, ^18O

    # ^*.
    char = re.sub('\^0' , '\u2070', char)
    char = re.sub('\^1' , '\u00B9', char)
    char = re.sub('\^2' , '\u00B2', char)
    char = re.sub('\^3' , '\u00B3', char)
    char = re.sub('\^4' , '\u2074', char)
    char = re.sub('\^5' , '\u2075', char)
    char = re.sub('\^6' , '\u2076', char)
    char = re.sub('\^7' , '\u2077', char)
    char = re.sub('\^8' , '\u2078', char)
    char = re.sub('\^9' , '\u2079', char)
    char = re.sub('\^0' , '\u2070', char)
    
    # ^**
    char = upper_upper_unicode(char)
    # ^***
    char = upper_upper_unicode(char)


    char = re.sub('1' , '\u2081', char)
    char = re.sub('2' , '\u2082', char)
    char = re.sub('3' , '\u2083', char)
    char = re.sub('4' , '\u2084', char)
    char = re.sub('5' , '\u2085', char)
    char = re.sub('6' , '\u2086', char)
    char = re.sub('7' , '\u2087', char)
    char = re.sub('8' , '\u2088', char)
    char = re.sub('9' , '\u2089', char)
    char = re.sub('\+' , '\u207A', char)
    char = re.sub(' \u207A ' , ' + ', char)
    char = re.sub('\-' , '\u207B', char)
    char = re.sub(' \u207B ' , ' - ', char)
    char = re.sub('\(\u2081D\)' , '(\u00B9D)', char)
    char = re.sub('\(\u2081S\)' , '(\u00B9S)', char)
    char = re.sub('\(\u2082D\)' , '(\u00B2D)', char)
    char = re.sub('\(\u2082P\)' , '(\u00B2P)', char)
    char = re.sub('\(\u2083P\)' , '(\u00B3P)', char)
    char = re.sub('\(\u2084P\)' , '(\u2074P)', char)
    char = re.sub('\(\u2084S\)' , '(\u2074S)', char)
    char = re.sub('\(v>\=\u2082\)' , '(\u03BD\u22652)', char)
    char = re.sub('\(v>\=\u2084\)' , '(\u03BD\u22654)', char)
    char = re.sub('hv' , 'h\u03BD', char)

    # _**
    char = lower_lower_unicode(char)
    # _***
    char = lower_lower_unicode(char)
    
    # ex) " ... nH2O"
    if nlabel0 >= 1:
        char = re.sub('\?na\?', num0[0], char)
    if nlabel0 >= 2:
        char = re.sub('\?nb\?', num0[1], char)
    if nlabel0 >= 3:
        char = re.sub('\?nc\?', num0[2], char)
    if nlabel0 >= 4:
        char = re.sub('\?nd\?', num0[3], char)
    if nlabel0 >= 5:
        char = re.sub('\?ne\?', num0[4], char)
    if nlabel0 >= 6:
        char = re.sub('\?nf\?', num0[5], char)
    if nlabel0 >= 7:
        char = re.sub('\?ng\?', num0[6], char)
    if nlabel0 >= 8:
        char = re.sub('\?nh\?', num0[7], char)
    if nlabel0 >= 9:
        char = re.sub('\?ni\?', num0[8], char)
    if nlabel0 >= 10:
        char = re.sub('\?nj\?', num0[9], char)

    # ex) "nH2O ...
    if nlabel1 == 1:
        char = re.sub('\?ma\?', num1[0], char)

    return char

# rate eq. -> unicode for display
def rate_unicode(rate):
    rate = re.sub('\*' , '\u00D7', rate)
    # 10^-**
    rate = re.sub('e-1' , ' \u00D7 10\u207B\u00B9', rate)
    rate = re.sub('e-2' , ' \u00D7 10\u207B\u00B2', rate)
    rate = re.sub('e-3' , ' \u00D7 10\u207B\u00B3', rate)
    rate = re.sub('e-4' , ' \u00D7 10\u207B\u2074', rate)
    rate = re.sub('e-5' , ' \u00D7 10\u207B\u2075', rate)
    rate = re.sub('e-6' , ' \u00D7 10\u207B\u2076', rate)
    rate = re.sub('e-7' , ' \u00D7 10\u207B\u2077', rate)
    rate = re.sub('e-8' , ' \u00D7 10\u207B\u2078', rate)
    rate = re.sub('e-9' , ' \u00D7 10\u207B\u2079', rate)
    # 10^**
    rate = re.sub('e1' , ' \u00D7 10\u00B9', rate)
    rate = re.sub('e2' , ' \u00D7 10\u00B2', rate)
    rate = re.sub('e3' , ' \u00D7 10\u00B3', rate)
    rate = re.sub('e4' , ' \u00D7 10\u2074', rate)
    rate = re.sub('e5' , ' \u00D7 10\u2075', rate)
    rate = re.sub('e6' , ' \u00D7 10\u2076', rate)
    rate = re.sub('e7' , ' \u00D7 10\u2077', rate)
    rate = re.sub('e8' , ' \u00D7 10\u2078', rate)
    rate = re.sub('e9' , ' \u00D7 10\u2079', rate)
    # ^-*.
    rate = re.sub('\^\(-0.' , '\u207B\u2070\u00B7', rate)
    rate = re.sub('\^\(-1.' , '\u207B\u00B9\u00B7', rate)
    rate = re.sub('\^\(-2.' , '\u207B\u00B2\u00B7', rate)
    rate = re.sub('\^\(-3.' , '\u207B\u00B3\u00B7', rate)
    rate = re.sub('\^\(-4.' , '\u207B\u2074\u00B7', rate)
    rate = re.sub('\^\(-5.' , '\u207B\u2075\u00B7', rate)
    rate = re.sub('\^\(-6.' , '\u207B\u2076\u00B7', rate)
    rate = re.sub('\^\(-7.' , '\u207B\u2077\u00B7', rate)
    rate = re.sub('\^\(-8.' , '\u207B\u2078\u00B7', rate)
    rate = re.sub('\^\(-9.' , '\u207B\u2079\u00B7', rate)
    rate = re.sub('\^\(-0.' , '\u207B\u2070\u00B7', rate)
    # ^*.
    rate = re.sub('\^0.' , '\u2070\u00B7', rate)
    rate = re.sub('\^1.' , '\u00B9\u00B7', rate)
    rate = re.sub('\^2.' , '\u00B2\u00B7', rate)
    rate = re.sub('\^3.' , '\u00B3\u00B7', rate)
    rate = re.sub('\^4.' , '\u2074\u00B7', rate)
    rate = re.sub('\^5.' , '\u2075\u00B7', rate)
    rate = re.sub('\^6.' , '\u2076\u00B7', rate)
    rate = re.sub('\^7.' , '\u2077\u00B7', rate)
    rate = re.sub('\^8.' , '\u2078\u00B7', rate)
    rate = re.sub('\^9.' , '\u2079\u00B7', rate)
    rate = re.sub('\^0.' , '\u2070\u00B7', rate)

    rate = re.sub('\^10.' , '\u00B9\u2070\u00B7', rate)
    rate = re.sub('\^11.' , '\u00B9\u00B9\u00B7', rate)
    rate = re.sub('\^12.' , '\u00B9\u00B2\u00B7', rate)
    rate = re.sub('\^13.' , '\u00B9\u00B3\u00B7', rate)
    rate = re.sub('\^14.' , '\u00B9\u2074\u00B7', rate)
    rate = re.sub('\^15.' , '\u00B9\u2075\u00B7', rate)
    rate = re.sub('\^16.' , '\u00B9\u2076\u00B7', rate)
    rate = re.sub('\^17.' , '\u00B9\u2077\u00B7', rate)
    rate = re.sub('\^18.' , '\u00B9\u2078\u00B7', rate)
    rate = re.sub('\^19.' , '\u00B9\u2079\u00B7', rate)

    # .*
    rate = re.sub('\u00B70' , '\u00B7\u2070', rate)
    rate = re.sub('\u00B71' , '\u00B7\u00B9', rate)
    rate = re.sub('\u00B72' , '\u00B7\u00B2', rate)
    rate = re.sub('\u00B73' , '\u00B7\u00B3', rate)
    rate = re.sub('\u00B74' , '\u00B7\u2074', rate)
    rate = re.sub('\u00B75' , '\u00B7\u2075', rate)
    rate = re.sub('\u00B76' , '\u00B7\u2076', rate)
    rate = re.sub('\u00B77' , '\u00B7\u2077', rate)
    rate = re.sub('\u00B78' , '\u00B7\u2078', rate)
    rate = re.sub('\u00B79' , '\u00B7\u2079', rate)

    # ^**
    rate = upper_upper_unicode(rate)
    # ^***
    rate = upper_upper_unicode(rate)

    rate = re.sub('\u00B9\)' , '\u00B9', rate)
    rate = re.sub('\u00B2\)' , '\u00B2', rate)
    rate = re.sub('\u00B3\)' , '\u00B3', rate)
    rate = re.sub('\u2074\)' , '\u2074', rate)
    rate = re.sub('\u2075\)' , '\u2075', rate)
    rate = re.sub('\u2076\)' , '\u2076', rate)
    rate = re.sub('\u2077\)' , '\u2077', rate)
    rate = re.sub('\u2078\)' , '\u2078', rate)
    rate = re.sub('\u2079\)' , '\u2079', rate)
    rate = re.sub('\u2070\)' , '\u2070', rate)

    rate = re.sub('k0' , 'k\u2080', rate)
    rate = re.sub('kinf' , 'k\u221E', rate)

    return rate

def unicode_LaTeX(char):

    # upper characters
    char = re.sub('\u00B9' , '$^1$', char)
    char = re.sub('\u00B2' , '$^2$', char)
    char = re.sub('\u00B3' , '$^3$', char)
    char = re.sub('\u2074' , '$^4$', char)
    char = re.sub('\u2075' , '$^5$', char)
    char = re.sub('\u2076' , '$^6$', char)
    char = re.sub('\u2077' , '$^7$', char)
    char = re.sub('\u2078' , '$^8$', char)
    char = re.sub('\u2079' , '$^9$', char)
    char = re.sub('\u2070' , '$^0$', char)

    # lower characters
    char = re.sub('\u2081' , '$_1$', char)
    char = re.sub('\u2082' , '$_2$', char)
    char = re.sub('\u2083' , '$_3$', char)
    char = re.sub('\u2084' , '$_4$', char)
    char = re.sub('\u2085' , '$_5$', char)
    char = re.sub('\u2086' , '$_6$', char)
    char = re.sub('\u2087' , '$_7$', char)
    char = re.sub('\u2088' , '$_8$', char)
    char = re.sub('\u2089' , '$_9$', char)
    char = re.sub('\u2080' , '$_0$', char)

    # upper dot
    char = re.sub('\u00B7' , '$^.$', char)

    # arrow
    char = re.sub('\u2192' , '$\\rightarrow$', char)

    # +, -, x
    char = re.sub('\u207A', '$^+$', char)
    char = re.sub('\u207B', '$^-$', char)
    char = re.sub('\u00D7', '$\\\\times$', char)

    # <, >
    char = re.sub('<', '$<$', char)
    char = re.sub('>', '$>$', char)

    # hv
    char = re.sub('hv', 'h$\\\\nu$', char)

    # Tn, Ti, Te
    char = re.sub('Tn', '$T_n$', char)
    char = re.sub('Ti', '$T_i$', char)
    char = re.sub('Te', '$T_e$', char)

    # k0, kinf
    char = re.sub('k0', '$k_0$', char)
    char = re.sub('kinf', '$k_\\\\infty$', char)

    return char

            
# Infix Notation to Reversed Poland Notation Parser #######################################################################
# split characters into tokens
def split_infix_tokens(eq):

    # insert ' ' among tokens
    #eq = eq.replace(' ','')
    eq = '('+eq+')'
    eq = eq.replace('T(electron)', 'Te')
    eq = eq.replace('T(ion)', 'Ti')
    eq = eq.replace('T(neutral)', 'Tn')
    eq = eq.replace('e-', 'eminus')
    eq = eq.replace('e+', 'eplus')
    eq = eq.replace('E-', 'eminus')
    eq = eq.replace('E+', 'eplus')
    eq = eq.replace('+', ' + ')
    eq = eq.replace('-', ' - ')
    eq = eq.replace('( + ' , '(+')
    eq = eq.replace('( - ' , '(-')
    eq = eq.replace('**', '^')
    eq = eq.replace('*', ' * ')
    eq = eq.replace('/', ' / ')
    eq = eq.replace('^', ' ^ ')
    eq = eq.replace('exp', ' exp ')
    eq = eq.replace('sqrt', ' sqrt ')
    eq = eq.replace('EXP', ' exp ')
    eq = eq.replace('SQRT', ' sqrt ')
    eq = eq.replace('(', ' ( ')
    eq = eq.replace(')', ' ) ')
    eq = eq.replace('eminus', 'e-')
    eq = eq.replace('eplus', 'e+')

    # split equation by ' '
    eq = eq.split(' ')

    out = []
    for i in range(len(eq)):
        eq[i] = eq[i].lstrip()
        eq[i] = eq[i].rstrip()
        if eq[i] != '':
            out.append(eq[i])
    return out

# translate infix into RPN
def infix_to_RPN(tokens):

    # constants for operator
    Left_assoc  = 0
    Right_assoc = 1

    # operators
    Operators = {
        '+'    : (0, Left_assoc),
        '-'    : (0, Left_assoc),
        '*'    : (1, Left_assoc),
        '/'    : (1, Left_assoc),
        '%'    : (1, Left_assoc),
        '^'    : (2, Right_assoc),
        'exp'  : (2, Right_assoc),
        'sqrt' : (2, Right_assoc)
    }

    out = []
    stack = []
    for token in tokens:
        if token in Operators.keys():
            while len(stack) != 0 and stack[-1] in Operators.keys():
                if (Operators[token][1] == Left_assoc  and Operators[token][0] - Operators[stack[-1]][0] <= 0) or(
                    Operators[token][1] == Right_assoc and Operators[token][0] - Operators[stack[-1]][0] <  0):
                    out.append(stack.pop())
                    continue
                break
            stack.append(token)
        elif token == '(':
            stack.append(token)
        elif token == ')':
            while len(stack) != 0 and stack[-1] != '(':
                out.append(stack.pop())
            stack.pop()
        else:
            out.append(token)
    while len(stack) != 0:
        out.append(stack.pop())
    return out

# for F90 input
def RPN_for_F90(token):

    # labels for the operators and terms
    Labels = {
        '+'    : (1, 1),
        '-'    : (1, 2),
        '*'    : (1, 3),
        '/'    : (1, 4),
        '%'    : (1, 5),
        '^'    : (1, 6),
        'exp'  : (1, 7),
        'sqrt' : (1, 8),
        'Tn'   : (2, 1),
        'Ti'   : (2, 2),
        'Te'   : (2, 3), 
        'M'    : (3, 1), 
        'h'    : (4, 1), 
    }

    out = [0]
    label = [0]
    for i in range(len(token)):
        if token[i] not in Labels.keys():
            out.append(token[i])
            label.append(0)
        else:
            out.append(Labels[token[i]][1])
            label.append(Labels[token[i]][0])
    out[0]   = len(out)-1
    label[0] = len(label)-1

    return out,label
            
# Analysis ################################################################################################################
#action: reaction analysis or output, action True or False?
#iplnt: which planet?
#reaction_chk_bln: checked reaction list of True or False
#fix_species_bln: select fixed species list
#input_species_char:
#fname_input_species:
def reaction_analysis(action, iplnt, reaction_chk_bln, fix_species_bln, dir0):

    print('Analyzing reactions...')

    Planet = Planet_list[iplnt][0]
    input_species_char = []
    input_species_path = []

    # Analysis of reaction_list
    reaction_list = []
    reaction_list_L = []
    reaction_list_R = []
    reaction_type_list = []
    reactant = []
    product = []
    rate_eq = []
    rate_eq_list = []
    reference_label = []
    reference = []
    label = []
    ref_label=[]
    species = []
    charge = []


    #reaction type
    #here is used to analyze reaction
    #normal reaction is set to 0 later
    type_photoionization  =  1 #reaction including hv in left side
    type_photodissociaion =  2 #reaction including hv in right side
    type_hv_right         =  4 #reaction including hv in right side
    type_e_impact         = 11 #electron impact reaction
    type_p_impact         = 12 #proton impact reaction
    type_H_impact         = 13 #H atom impact reaction
    type_pressure_3body   = 30 #3 body reaction
    type_pressure_3bodyM  = 31 #3 body reaction: k*M
    type_Lindemann_Hinshelwood = 32
    type_Meteoroid        = 41 #Meteoroid ablation
    type_Rainout          = 47 #Rainout

    #this index counts selected reaction number in the loop below
    isp = 0
    ich = 0
    #loop reaction list
    for i in range(len(reaction_rate_list)):
        if reaction_chk_bln[i].get() == True: #only checked reaction

            #print(reaction_rate_list[i]) # for confirmation

            # treatment of "X+ + X- -> products" reaction
            # replace it with ion-ion recombination reacrtions "ion+ + X- -> products" and "ion- + X+ -> products"


            # Analysis of reactant and product
            #reaction->  "reaction : coefficient"
            #findall function returns as list!
            reaction_left_right = re.findall('(.*):', reaction_rate_list[i])  #finding : and get "reaction(string type) part" ,removing coefficient part

            reaction_list_tmp = reaction_left_right[0]
            reaction_list_tmp = reaction_unicode(reaction_list_tmp)
            reaction_list.append(reaction_list_tmp)

            #left -> right reaction
            #separating reaction into left(reactants) and right(products)
            reaction_left = re.findall('(.*)->', reaction_left_right[0])
            reaction_list_L.append(reaction_left[0])
            reaction_right = re.findall('->(.*)', reaction_left_right[0])
            reaction_list_R.append(reaction_right[0])

            itype = 0

            # type determination
            if 'products' in reaction_right[0]:
                reaction_right[0] = reaction_right[0].strip(' products ')
                reaction_right[0] = reaction_right[0].strip('+')
            if 'Condensation' in reaction_right[0]:
                reaction_right[0] = reaction_right[0].strip(' Condensation ')
                reaction_right[0] = reaction_right[0].strip('+')
            # photoionization, photodisosciation, excitation
            if 'hv' in reaction_left[0]:
                if 'Photoionization' in reaction_rate_list[i]:
                    itype = type_photoionization
                    reaction_left[0] = reaction_left[0].strip(' hv ')
                    reaction_left[0] = reaction_left[0].strip('+')
                if 'Photodissociation' in reaction_rate_list[i]:
                    itype = type_photodissociaion
                    reaction_left[0] = reaction_left[0].strip(' hv ')
                    reaction_left[0] = reaction_left[0].strip('+')
            # airglow
            if 'hv' in reaction_right[0]:
                itype = type_hv_right
                reaction_right[0] = reaction_right[0].strip(' hv ')
                reaction_right[0] = reaction_right[0].strip('+')
            # pressure dependent three body reactions
            if 'M' in reaction_left[0] and 'k0' in reaction_rate_list[i]:
                itype = type_pressure_3body
                if 'M' in reaction_right[0]:
                    reaction_right[0] = reaction_right[0].strip(' M ')
                    reaction_right[0] = reaction_right[0].strip('+')
            # pressure dependent three body reactions
            if 'M' in reaction_left[0] and 'k0M' in reaction_rate_list[i]:
                itype = type_pressure_3bodyM
                if 'M' in reaction_right[0]:
                    reaction_right[0] = reaction_right[0].strip(' M ')
                    reaction_right[0] = reaction_right[0].strip('+')
            # Lindemann-Hinshelwood expression
            if 'k0LH' in reaction_rate_list[i]:
                itype = type_Lindemann_Hinshelwood
            # electron impact reactions
            if 'e-*' in reaction_left[0]:
                itype = type_e_impact
                reaction_left[0]= reaction_left[0].strip(' e-* ')
                reaction_left[0] = reaction_left[0].strip('+')
                reaction_right[0] = reaction_right[0].strip(' e-* ')
                reaction_right[0] = reaction_right[0].rstrip('+')
            # proton impact reactions
            if 'p*' in reaction_left[0]:
                itype = type_p_impact
                reaction_left[0] = reaction_left[0].strip(' p* ')
                reaction_left[0] = reaction_left[0].strip('+')
                reaction_right[0] = reaction_right[0].strip(' p* ')
                reaction_right[0] = reaction_right[0].rstrip('+')
            # H impact reactions
            if 'H*' in reaction_left[0]:
                itype = type_p_impact
                reaction_left[0] = reaction_left[0].strip(' H* ')
                reaction_left[0] = reaction_left[0].strip('+')
                reaction_right[0] = reaction_right[0].strip(' H* ')
                reaction_right[0] = reaction_right[0].rstrip('+')
            # Meteoroid ablation
            if 'Meteoroid' in reaction_left[0]:
                itype = type_Meteoroid
                reaction_left[0] = reaction_left[0].strip(' Meteoroid ')
                reaction_left[0] = reaction_left[0].strip('+')
            if 'Rainout' in reaction_right[0]:
                itype = type_Rainout
                reaction_right[0] = reaction_right[0].strip(' Rainout ')
                reaction_right[0] = reaction_right[0].strip('+')


            #reactant = list of reaction_left
            #reactant
            #[['sp1','sp2',,,], [], [], [], [] ,,,]
            reactant.append(reaction_left[0].split(' + '))

            #to remove all spaces
            #ex) " CO2  " -> "CO2"
            for j in range(len(reactant[ich])):
                reactant[ich][j] = reactant[ich][j].lstrip()
                reactant[ich][j] = reactant[ich][j].rstrip()
                if re.compile('^\d+?').search(reactant[ich][j]):
                    num = re.findall('^(\d{1,3})', reactant[ich][j])
                    reactant[ich][j] = re.sub('^\d{1,3}', '', reactant[ich][j])
                    for k in range(int(num[0])-1):
                        reactant[ich].append(reactant[ich][j])

            product.append(reaction_right[0].split(' + '))
            for j in range(len(product[ich])):
                product[ich][j] = product[ich][j].lstrip()
                product[ich][j] = product[ich][j].rstrip()
                if re.compile('^\d+?').search(product[ich][j]):
                    num = re.findall('^(\d{1,3})', product[ich][j])
                    product[ich][j] = re.sub('^\d{1,3}', '', product[ich][j])
                    for k in range(int(num[0])-1):
                        product[ich].append(product[ich][j])
            #print(product[ich])

            reaction_type_list.append(itype)

            # Analysis of Species
            for j in range(len(reactant[ich])):
                label = 0
                for isp in range(len(species)):
                    if reactant[ich][j] == species[isp]:
                        label = 1
                if label == 0 and reactant[ich][j] != '':
                    species.append(reactant[ich][j])

            for j in range(len(product[ich])):
                label = 0
                for isp in range(len(species)):
                    if product[ich][j] == species[isp]:
                        label = 1
                if label == 0 and product[ich][j] != '':
                    species.append(product[ich][j])

            # Analysis of reaction rate (not detailed) : detailed analysis > later
            # getting rate coefficient part, removing reaction part in the left
            # this returns as a list type!!
            rate_ref_label_tmp = re.findall(':(.*)', reaction_rate_list[i])

            # @ denotes reference and # denotes label. Some reaction have reference with coefficient
            rate_ref_label_all_cases = rate_ref_label_tmp[0].split('&&')
            rate_all_cases = ['' for i in range(len(rate_ref_label_all_cases))]
            ref_all_cases = ['' for i in range(len(rate_ref_label_all_cases))]
            label_all_cases = ['' for i in range(len(rate_ref_label_all_cases))]
            for i in range(len(rate_ref_label_all_cases)):
                if '@' not in rate_ref_label_all_cases[i] and '#' not in rate_ref_label_all_cases[i]:
                    rate_all_cases[i] = rate_ref_label_all_cases[i]
                if '@' in rate_ref_label_all_cases[i] and '#' not in rate_ref_label_all_cases[i]:
                    rate_tmp = re.findall('(.*)@', rate_ref_label_all_cases[i])
                    rate_all_cases[i] = rate_tmp[0]
                    ref_tmp  = re.findall('@(.*)', rate_ref_label_all_cases[i])
                    ref_all_cases[i] = ref_tmp[0]
                if '@' not in rate_ref_label_all_cases[i] and '#' in rate_ref_label_all_cases[i]:
                    rate_tmp = re.findall('(.*)#', rate_ref_label_all_cases[i])
                    rate_all_cases[i] = rate_tmp[0]
                    label_tmp  = re.findall('#(.*)', rate_ref_label_all_cases[i])
                    label_all_cases[i] = label_tmp[0]
                if '@' in rate_ref_label_all_cases[i] and '#' in rate_ref_label_all_cases[i]:
                    rate_tmp = re.findall('(.*)@', rate_ref_label_all_cases[i])
                    if '#' in rate_tmp[0]:
                        rate_tmp1 = re.findall('(.*)#', rate_tmp[0])
                        rate_all_cases[i] = rate_tmp1[0]
                        ref_tmp  = re.findall('@(.*)', rate_ref_label_all_cases[i])
                        ref_all_cases[i] = ref_tmp[0]
                        label_tmp  = re.findall('#(.*)@', rate_ref_label_all_cases[i])
                        label_all_cases[i] = label_tmp[0]
                    if '#' not in rate_tmp[0]:
                        rate_all_cases[i] = rate_tmp[0]
                        ref_tmp  = re.findall('@(.*)#', rate_ref_label_all_cases[i])
                        ref_all_cases[i] = ref_tmp[0]
                        label_tmp  = re.findall('#(.*)', rate_ref_label_all_cases[i])
                        label_all_cases[i] = label_tmp[0]
            rate_eq.append(list(rate_all_cases))
            rate_eq_list.append(list(rate_all_cases))
            reference.append(list(ref_all_cases))
            ref_label.append(list(label_all_cases))

            # update checked index
            ich = ich + 1

    #to check all the reaction have type
    if len(reaction_type_list) != len(reactant):
        print('reaction type error', len(reaction_type_list), len(reactant))

    # reactant, product list
    #[[0,0,0,0,0,0,0,0,0,0],[],[],,,]
    reactant_list = [ [0 for isp in range(100)] for ich in range(len(reactant))]
    product_list  = [ [0 for isp in range(100)] for ich in range(len(product))]


    for ich in range(len(reactant)):
        reactant_list[ich][0] = len(reactant[ich])
        for isp in range(len(reactant[ich])):
            for jsp in range(len(species)):
                if species[jsp] == reactant[ich][isp]:
                    reactant_list[ich][isp+1] = jsp+1
        product_list[ich][0] = len(product[ich])
        for isp in range(len(product[ich])):
            for jsp in range(len(species)):
                if species[jsp] == product[ich][isp]:
                    product_list[ich][isp+1] = jsp+1

    # Check for the reactions identified as the same
    same_label = 0
    txt = ''
    for ich in range(len(reactant)):

        type_i = reaction_type_list[ich]

        rlist_i = [0 for i in range(len(reactant[ich]))]
        for i in range(len(reactant[ich])):
            rlist_i[i] = reactant_list[ich][i+1]
        rlist_i.sort()

        plist_i = [0 for i in range(len(product[ich]))]
        for i in range(len(product[ich])):
            plist_i[i] = product_list[ich][i+1]
        plist_i.sort()

        for jch in range(len(reactant)-ich-1):

            type_j = reaction_type_list[jch+ich+1]

            rlist_j = [0 for j in range(len(reactant[jch+ich+1]))]
            for j in range(len(reactant[jch+ich+1])):
                rlist_j[j] = reactant_list[jch+ich+1][j+1]
            rlist_j.sort()

            plist_j = [0 for j in range(len(product[jch+ich+1]))]
            for j in range(len(product[jch+ich+1])):
                plist_j[j] = product_list[jch+ich+1][j+1]
            plist_j.sort()

            if rlist_i == rlist_j and plist_i == plist_j and type_i == type_j:
                same_label = 1
                txt+='error! Following reactions are identified as the same reaction.\n'
                txt+='  '+ reaction_list[ich]
                for i in range(len(rate_eq[ich])):
                    if reference[ich][i] != '' and ref_label[ich][i] != '':
                        txt+='     '+rate_eq[ich][i]+' @ '+reference[ich][i]+' # '+ref_label[ich][i]+'\n'
                    if reference[ich][i] == '' and ref_label[ich][i] != '':
                        txt+='     '+rate_eq[ich][i]+' # '+ref_label[ich][i]+'\n'
                    if reference[ich][i] != '' and ref_label[ich][i] == '':
                        txt+='     '+rate_eq[ich][i]+' @ '+reference[ich][i]+'\n'
                    if reference[ich][i] == '' and ref_label[ich][i] == '':
                        txt+='     '+rate_eq[ich][i]+'\n'
                txt+='  '+ reaction_list[jch+ich+1]
                for i in range(len(rate_eq[jch+ich+1])):
                    if reference[jch+ich+1][i] != '' and ref_label[jch+ich+1][i] != '':
                        txt+='     '+rate_eq[jch+ich+1][i]+' @ '+reference[jch+ich+1][i]+' # '+ref_label[jch+ich+1][i]+'\n'
                    if reference[jch+ich+1][i] == '' and ref_label[jch+ich+1][i] != '':
                        txt+='     '+rate_eq[jch+ich+1][i]+' # '+ref_label[jch+ich+1][i]+'\n'
                    if reference[jch+ich+1][i] != '' and ref_label[jch+ich+1][i] == '':
                        txt+='     '+rate_eq[jch+ich+1][i]+' @ '+reference[jch+ich+1][i]+'\n'
                    if reference[jch+ich+1][i] == '' and ref_label[jch+ich+1][i] == '':
                        txt+='     '+rate_eq[jch+ich+1][i]+'\n'
                txt+='\n'

    if same_label == 1:
        action = 'error'
        err_win = tk.Toplevel() #display on the main window
        err_win.title("error")
        err_win.geometry("720x640")
        txt_label = tk.Label(err_win, anchor="w", justify="left", font=('', 15), text=txt)
        txt_label.place(x=10, y=10)


    # input species
    #0 -> not input, so input 0 in this case
    #1 -> input
    # read path of fixed species
    input_species = [0] * len(species)
    fixbln = [0] * len(species)
    path = './'+Planet+'/'+dir0+'/settings/input_path_list.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            line = f.readline()
            input_sp = line.split(',')
            input_dir = f.readline().strip('\n')
            sp_fname = f.readlines()
            for isp in range(len(input_sp)):
                input_sp[isp] = input_sp[isp].lstrip().rstrip()
                for jsp in range(len(sp_fname)):
                    sp = re.findall('(.*):',sp_fname[jsp])
                    sp = sp[0].lstrip().rstrip()
                    if input_sp[isp] == sp:
                        input_species_char.append(sp)
                        fname = re.findall(':(.*)',sp_fname[jsp])
                        fname = fname[0].strip('\n').lstrip().rstrip()
                        if ';' in fname:
                            fixbln[isp] = re.findall(';(.*)', fname)[0].lstrip().rstrip()
                            fname = re.findall('(.*);', fname)[0]
                        input_species_path.append('./'+Planet+'/'+input_dir+'/'+fname)
        f.close()
        for isp in range(len(input_species_char)):
            input_species_path[isp].replace('//','/')
            for jsp in range(len(species)):
                if species[jsp] == sp:
                    input_species[jsp] = 1

    if os.path.exists(path) == False:
        err_win = tk.Toplevel() #display on the main window
        err_win.title("error")
        err_win.geometry("300x100")
        txt= path+'" does not exist!\nPlease set input species and their paths again!'
        txt_label = tk.Label(err_win, anchor="w", justify="left", font=('', 15), text=txt)
        txt_label.place(x=10, y=10)

    # fixed species
    #0 -> variable
    #1 -> fixed
    fix_species = [0] * len(species)
    for isp in range(len(species)):
        for jsp in range(len(input_species_char)):
            if species[isp] == input_species_char[jsp]:
                if fixbln[jsp] == '1':
                    fix_species[isp] = 1
        if species[isp] == 'M':
            fix_species[isp] = 1
        losslabel = 0
        for ich in range(len(reactant)):
            for j in range(len(reactant[ich])):
                if species[isp] == reactant[ich][j]:
                    losslabel = 1
        # If a species has no loss reactions: it is fixed
        if losslabel == 0:
            fix_species[isp] = 1

    # boundary condition
    # make type of boundary condition check button
    bc_read = {}
    bc_species = {}
    bc_Llabel  = {}
    bc_Lval    = {}
    bc_Ulabel  = {}
    bc_Uval    = {}
    # read path of fixed species
    path = './'+Planet+'/'+dir0+'/settings/boundary_condition.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            line = f.readline()
            bc_read = f.readlines()
            for isp in range(len(bc_read)):
                bc_read[isp] = bc_read[isp].strip('\n')
                bc_species[isp] = re.findall('(.*):',bc_read[isp])[0].lstrip().rstrip()
                bcchar = re.findall(':(.*)',bc_read[isp])
                bc_Llabel[isp] = re.findall('L(.*)lval',bcchar[0])[0].lstrip().rstrip()
                bc_Lval[isp]   = re.findall('lval(.*)U',bcchar[0])[0].lstrip().rstrip()
                bc_Ulabel[isp] = re.findall('U(.*)uval',bcchar[0])[0].lstrip().rstrip()
                bc_Uval[isp]   = re.findall('uval(.*)',bcchar[0])[0].lstrip().rstrip()
        f.close()

    # Charge and Mass definition
    #
    mass = [0.00000] * len(species)
    for isp in range(len(species)):

        # Charge
        # calc charge of each species
        charge.append(str(species[isp].count("+") - species[isp].count("-")))

        ##### Mass calculation ################################################
        # remove unnecessary symbols
        species_tmp = species[isp]
        species_tmp = re.sub('\(1D\)', '', species_tmp)
        species_tmp = re.sub('\(1S\)', '', species_tmp)
        species_tmp = re.sub('\(2D\)', '', species_tmp)
        species_tmp = re.sub('\(2P\)', '', species_tmp)
        species_tmp = re.sub('\(3P\)', '', species_tmp)
        species_tmp = re.sub('\(4P\)', '', species_tmp)
        species_tmp = re.sub('\(4S\)', '', species_tmp)
        species_tmp = re.sub('\(v>=2\)', '', species_tmp)
        species_tmp = re.sub('\(v>=4\)', '', species_tmp)
        species_tmp = re.sub('\+', '', species_tmp) #string until here
        species_tmp = re.sub('\-', '', species_tmp) #string until here

        # ex) (H2O)2 -> H2OH2O
        if re.compile('\(.*\)\d+?').search(species_tmp): #find spcies including (***)n
            species_brackets = re.findall('\((.*)\)\d+?', species_tmp) #find spcies in brackets
            num_str = re.findall('\)(\d+?)', species_tmp) #find number outside the brackets
            num_int = int(num_str[0], base = 10) #string -> number, base=10: demical number
            for k in range(num_int-1): #loop "number" times, ex) (H2O)2 -> 2 times, (H2O)3 -> 3 times
                species_tmp = re.sub('\)(\d+?)', '', species_tmp) #delete number, ex (H2O)2-> O
                species_tmp = re.sub('\(', '', species_tmp) #delete (
                species_tmp = re.sub('\)', '', species_tmp) #delete )
                species_tmp = species_tmp + species_brackets[0] #ex) append H2O
            #print(species_tmp)

        nisotope = 0
        isotopes = []
        if re.compile('\^\d+?').search(species_tmp):
            isotopes = re.findall('\^\d{1,3}[A-Z][^A-Z|\^]*', species_tmp)
            nisotope = len(isotopes)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-aa', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ab', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ac', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ad', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ae', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-af', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ag', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ah', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ai', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-aj', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ak', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-al', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-am', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-an', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ao', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ap', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-aq', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ar', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-as', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-at', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-au', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-av', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-aw', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ax', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ay', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-az', species_tmp,1)

            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-ba', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bb', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bc', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bd', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-be', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bf', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bg', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bh', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bi', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bj', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bk', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bl', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bm', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bn', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bo', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bp', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bq', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-br', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bs', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bt', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bu', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bv', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bw', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bx', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-by', species_tmp,1)
            species_tmp = re.sub('\^\d{1,3}[A-Z][^A-Z|\^]*', 'Isotope-bz', species_tmp,1)
            #print(species_tmp, isotopes)

        species_split = re.findall('[A-Z][^A-Z]*', species_tmp) #split elements

        if nisotope >= 1:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-aa', isotopes[0], species_split[i],1)
        if nisotope >= 2:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ab', isotopes[1], species_split[i],1)
        if nisotope >= 3:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ac', isotopes[2], species_split[i],1)
        if nisotope >= 4:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ad', isotopes[3], species_split[i],1)
        if nisotope >= 5:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ae', isotopes[4], species_split[i],1)
        if nisotope >= 6:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-af', isotopes[5], species_split[i],1)
        if nisotope >= 7:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ag', isotopes[6], species_split[i],1)
        if nisotope >= 8:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ah', isotopes[7], species_split[i],1)
        if nisotope >= 9:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ai', isotopes[8], species_split[i],1)
        if nisotope >= 10:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-aj', isotopes[9], species_split[i],1)
        if nisotope >= 11:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ak', isotopes[10], species_split[i],1)
        if nisotope >= 12:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-al', isotopes[11], species_split[i],1)
        if nisotope >= 13:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-am', isotopes[12], species_split[i],1)
        if nisotope >= 14:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-an', isotopes[13], species_split[i],1)
        if nisotope >= 15:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ao', isotopes[14], species_split[i],1)
        if nisotope >= 16:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ap', isotopes[15], species_split[i],1)
        if nisotope >= 17:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-aq', isotopes[16], species_split[i],1)
        if nisotope >= 18:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ar', isotopes[17], species_split[i],1)
        if nisotope >= 19:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-as', isotopes[18], species_split[i],1)
        if nisotope >= 20:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-at', isotopes[19], species_split[i],1)
        if nisotope >= 21:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-au', isotopes[20], species_split[i],1)
        if nisotope >= 22:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-av', isotopes[21], species_split[i],1)
        if nisotope >= 23:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-aw', isotopes[22], species_split[i],1)
        if nisotope >= 24:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ax', isotopes[23], species_split[i],1)
        if nisotope >= 25:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ay', isotopes[24], species_split[i],1)
        if nisotope >= 26:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-az', isotopes[25], species_split[i],1)

        if nisotope >= 27:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-ba', isotopes[26], species_split[i],1)
        if nisotope >= 28:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bb', isotopes[27], species_split[i],1)
        if nisotope >= 29:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bc', isotopes[28], species_split[i],1)
        if nisotope >= 30:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bd', isotopes[29], species_split[i],1)
        if nisotope >= 31:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-be', isotopes[30], species_split[i],1)
        if nisotope >= 32:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bf', isotopes[31], species_split[i],1)
        if nisotope >= 33:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bg', isotopes[32], species_split[i],1)
        if nisotope >= 34:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bh', isotopes[33], species_split[i],1)
        if nisotope >= 35:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bi', isotopes[34], species_split[i],1)
        if nisotope >= 36:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bj', isotopes[35], species_split[i],1)
        if nisotope >= 37:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bk', isotopes[36], species_split[i],1)
        if nisotope >= 38:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bl', isotopes[37], species_split[i],1)
        if nisotope >= 39:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bm', isotopes[38], species_split[i],1)
        if nisotope >= 40:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bn', isotopes[39], species_split[i],1)
        if nisotope >= 41:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bo', isotopes[40], species_split[i],1)
        if nisotope >= 42:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bp', isotopes[41], species_split[i],1)
        if nisotope >= 43:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bq', isotopes[42], species_split[i],1)
        if nisotope >= 44:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-br', isotopes[43], species_split[i],1)
        if nisotope >= 45:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bs', isotopes[44], species_split[i],1)
        if nisotope >= 46:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bt', isotopes[45], species_split[i],1)
        if nisotope >= 47:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bu', isotopes[46], species_split[i],1)
        if nisotope >= 48:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bv', isotopes[47], species_split[i],1)
        if nisotope >= 49:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bw', isotopes[48], species_split[i],1)
        if nisotope >= 50:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bx', isotopes[49], species_split[i],1)
        if nisotope >= 51:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-by', isotopes[50], species_split[i],1)
        if nisotope >= 52:
            for i in range(len(species_split)):
                species_split[i] = re.sub('Isotope-bz', isotopes[51], species_split[i],1)
        #print(species_split)

        for j in range(len(species_split)):
            num_str = re.findall('[A-Z]*(\d*)$', species_split[j]) #find spcies which start with capital and end with number
            # ex) CO2 -> C,O,O
            if num_str[0] != '':
                num_int = int(num_str[0], base = 10) #string -> number
                for k in range(num_int-1): #loop "number" times, ex) O2 -> 2 times, O3 -> 3 times
                    species_split[j] = re.sub('\d*$', '', species_split[j]) #delete number, ex O2-> O
                    species_split.append(species_split[j]) #ex) append O

        #print(species_split)

        # atom and molecule mass
        for j in range(len(species_split)):
            if species_split[j] == 'H':
                mass[isp] = mass[isp] + 1.0
            if species_split[j] == 'He':
                mass[isp] = mass[isp] + 4.0
            if species_split[j] == 'C':
                mass[isp] = mass[isp] + 12.0
            if species_split[j] == '^1C': # singlet C   e.g.) ^1CH3 
                mass[isp] = mass[isp] + 12.0
            if species_split[j] == '^3C': # triplet C   e.g.) ^3CH3 
                mass[isp] = mass[isp] + 12.0
            if species_split[j] == '^13C':
                mass[isp] = mass[isp] + 13.0
            if species_split[j] == '^14C':
                mass[isp] = mass[isp] + 14.0
            if species_split[j] == 'N':
                mass[isp] = mass[isp] + 14.0
            if species_split[j] == '^15N':
                mass[isp] = mass[isp] + 15.0
            if species_split[j] == 'O':
                mass[isp] = mass[isp] + 16.0
            if species_split[j] == '^17O':
                mass[isp] = mass[isp] + 17.0
            if species_split[j] == '^18O':
                mass[isp] = mass[isp] + 18.0
            if species_split[j] == 'Na':
                mass[isp] = mass[isp] + 23.0
            if species_split[j] == 'Mg':
                mass[isp] = mass[isp] + 24.0
            if species_split[j] == 'Si':
                mass[isp] = mass[isp] + 28.0
            if species_split[j] == 'S':
                mass[isp] = mass[isp] + 32.0
            if species_split[j] == 'Cl':
                mass[isp] = mass[isp] + 35.0
            if species_split[j] == 'Ar':
                mass[isp] = mass[isp] + 40.0
            if species_split[j] == 'Fe':
                mass[isp] = mass[isp] + 56.0
        # ion and electron mass
        mass[isp] = mass[isp] - 0.00054858 * float(species[isp].count("+") - species[isp].count("-"))
        # followings are fixed species, but avoiding 0 mass
        if species[isp] == 'M':
            mass[isp] = 10.00000
        if species[isp] == 'hv':
            mass[isp] = 10.00000
        if species[isp] == 'products':
            mass[isp] = 10.00000

    # Production, Loss and Jacobian Analysis ##################################
    # Production list
    #[ [number of reaction associated with production of sp,reaction index (i), the number of sp produced in (i),reaction index (j),the number of sp produced in (j),,,]  for sp1
    #  [,,,] for sp2
    #  [,,,] for sp3
    #  ,,,
    #  ,,,
    #  [,,,]   ]
    Production_label = [[0 for i in range(len(product))] for j in range(len(species))]
    Production_list = [[0] for j in range(len(species))]
    for isp in range(len(species)):
        nch = 0
        for ich in range(len(product)):
            for jsp in range(len(product[ich])):
                if species[isp] == product[ich][jsp]:
                    Production_label[isp][ich] = Production_label[isp][ich] + 1
                    for i in range(10):
                        if Production_label[isp][ich] == i+2:
                            Production_list[isp].remove([ich+1,i+1])
                            nch = nch - 1
                    Production_list[isp].append([ich+1,Production_label[isp][ich]])
                    nch = nch + 1
        Production_list[isp][0] = [nch]

    #search max
    # get maximum of Production_list[:][0][0]
    # this is needed to make an output for Fortran!
    nch_arr = [0]
    for isp in range(len(species)):
        if fix_species[isp] == 0:
            nch_arr.append(Production_list[isp][0][0])
    max_ch_P = max(nch_arr) * 2 + 1

    Production_list_for_output = [[-1 for i in range(max_ch_P)] for j in range(len(species))]
    for isp in range(len(species)):
        if fix_species[isp] == 0:
            for i in range(len(Production_list[isp])):
                if i == 0:
                    Production_list_for_output[isp][0] = Production_list[isp][0][0]
                if i != 0:
                    for j in range(len(Production_list[isp][i])):
                        Production_list_for_output[isp][2*i+j-1] = Production_list[isp][i][j]

    # Loss list
    #[ [number of reaction associated with loss of sp,reaction index (i), the number of sp lost in (i),reaction index (j),the number of sp lost in (j),,,]  for sp1
    #  [,,,] for sp2
    #  [,,,] for sp3
    #  ,,,
    #  ,,,
    #  [,,,]   ]
    Loss_label = [[0 for i in range(len(reactant))] for j in range(len(species))]
    Loss_list = [[0] for j in range(len(species))]
    for isp in range(len(species)):
        nch = 0
        for ich in range(len(reactant)):
            for jsp in range(len(reactant[ich])):
                if species[isp] == reactant[ich][jsp]:
                    Loss_label[isp][ich] = Loss_label[isp][ich] + 1
                    for i in range(10):
                        if Loss_label[isp][ich] == i+2:
                            Loss_list[isp].remove([ich+1,i+1])
                            nch = nch - 1
                    Loss_list[isp].append([ich+1,Loss_label[isp][ich]])
                    nch = nch + 1
        Loss_list[isp][0] = [nch]

    # search max
    # get maximum of Loss_list[:][0][0]
    # this is needed to make an output for Fortran!
    nch_arr = [0]
    for isp in range(len(species)):
        if fix_species[isp] == 0:
            nch_arr.append(Loss_list[isp][0][0])
    max_ch_L = max(nch_arr) * 2 + 1

    Loss_list_for_output = [[-1 for i in range(max_ch_L)] for j in range(len(species))]
    for isp in range(len(species)):
        if fix_species[isp] == 0:
            for i in range(len(Loss_list[isp])):
                if i == 0:
                    Loss_list_for_output[isp][0] = Loss_list[isp][0][0]
                if i != 0:
                    for j in range(len(Loss_list[isp][i])):
                        Loss_list_for_output[isp][2*i+j-1] = Loss_list[isp][i][j]

    #### Jacobian list #########################################################
    # Jacobian_list = [i,j, reaction index k,contribution of k, reaction index l, contribution of l,,,,,]
    # matrix (i,j)
    # derivative species i of species j
    # d species_i / d species_j

    #Jacobian_label = i*j matrix
    # 1-> has some variable
    # 0-> does not have variable
    Jacobian_label_P = [[0 for i in range(len(species))] for j in range(len(species))]
    Jacobian_label_L = [[0 for i in range(len(species))] for j in range(len(species))]
    Jacobian_list = [[[0] for j in range(len(species))] for k in range(len(species))]
    for isp in range(len(species)):
        for jsp in range(len(species)):
            nch = 0
            for ich in range(len(reactant)): #the number of reactions
                if Production_label[isp][ich] != 0: #if ith species have some variable
                    for ksp in range(len(reactant[ich])): #loop over reactant of ith species
                        if species[jsp] == reactant[ich][ksp]: #if jth species == ith species
                            Jacobian_label_P[isp][jsp] = 1
                            for i in range(len(Jacobian_list[isp][jsp])):
                                if Jacobian_list[isp][jsp][i] == [ich+1,Production_label[isp][ich]]:
                                    Jacobian_list[isp][jsp].remove([ich+1,Production_label[isp][ich]])
                                    nch = nch - 1
                            Jacobian_list[isp][jsp].append([ich+1,Production_label[isp][ich]])
                            nch = nch + 1
            for ich in range(len(reactant)):
                if Loss_label[isp][ich] != 0:
                    for ksp in range(len(reactant[ich])):
                        if species[jsp] == reactant[ich][ksp]:
                            Jacobian_label_L[isp][jsp] = 1
                            for i in range(len(Jacobian_list[isp][jsp])):
                                if Jacobian_list[isp][jsp][i] == [ich+1,-Loss_label[isp][ich]]:
                                    Jacobian_list[isp][jsp].remove([ich+1,-Loss_label[isp][ich]])
                                    nch = nch - 1
                            Jacobian_list[isp][jsp].append([ich+1,-Loss_label[isp][ich]])
                            nch = nch + 1
            Jacobian_list[isp][jsp][0] = [nch]

    nch_arr = [0]
    for isp in range(len(species)):
        if fix_species[isp] == 0:
            for jsp in range(len(species)):
                if fix_species[jsp] == 0:
                    nch_arr.append(Jacobian_list[isp][jsp][0][0])
    max_ch_J = max(nch_arr) * 2 + 1

    # i*j Jacobian matrix for output ######
    Jacobian_list_for_output = [[[0 for i in range(max_ch_J)] for j in range(len(species))] for k in range(len(species))]
    for isp in range(len(species)):
        if fix_species[isp] == 0:
            for jsp in range(len(species)):
                if fix_species[jsp] == 0:
                    for i in range(len(Jacobian_list[isp][jsp])):
                        if i == 0:
                            Jacobian_list_for_output[isp][jsp][0] = Jacobian_list[isp][jsp][0][0]
                        if i != 0:
                            for j in range(len(Jacobian_list[isp][jsp][i])):
                                Jacobian_list_for_output[isp][jsp][2*i+j-1] = Jacobian_list[isp][jsp][i][j]


    #### Analysis of reaction rate expression ##################################

    # Max number of cases of T dependence, three body reaction
    ncases = 3

    #rate_
    rate_cases  = [ 0 for ich in range(len(rate_eq))]
    rate_Trange = [[[0 for i in range(3)] for ieq in range(ncases)] for ich in range(len(rate_eq))]
    rate_rpn_token = [[[] for ieq in range(ncases)] for ich in range(len(rate_eq))]
    rate_rpn_label = [[[] for ieq in range(ncases)] for ich in range(len(rate_eq))]

    ntokens = 0
    for ich in range(len(rate_eq)):

        ################# T range array ######################
        # rate_Trange[ich][0] = [type, start [K], end[K]]
        # T range default
        if len(rate_eq[ich]) == 1:
            rate_Trange[ich][0][1] = 0
            rate_Trange[ich][0][2] = 100000

        #excluding p-dependent
        if reaction_type_list[ich] != type_pressure_3body or reaction_type_list[ich] != type_pressure_3bodyM or reaction_type_list[ich] != type_Lindemann_Hinshelwood:

            if len(rate_eq[ich]) >= 1:
                #T type
                #rate_eq[ich][ieq][0]: type of Temperature
                for ieq in range(len(rate_eq[ich])):
                    if 'Tn' in rate_eq[ich][ieq]:
                        rate_Trange[ich][ieq][0] = 1
                    if 'Te' in rate_eq[ich][ieq]:
                        rate_Trange[ich][ieq][0] = 2
                    if 'Ti' in rate_eq[ich][ieq]:
                        rate_Trange[ich][ieq][0] = 3
                    if 'T' in rate_eq[ich][ieq] and '~' in rate_eq[ich][ieq] and '[K]' in rate_eq[ich][ieq] and 'Tn' not in rate_eq[ich][ieq] and 'Te' not in rate_eq[ich][ieq] and 'Ti' not in rate_eq[ich][ieq]:
                        rate_Trange[ich][ieq][0] = 4

                #T type exclusion -> 4
                for ieq in range(len(rate_eq[ich])):
                    if rate_Trange[ich][ieq][0] == 4:
                        for jeq in range(len(rate_eq[ich])):
                            if rate_Trange[ich][jeq][0] != 4:
                                rate_Trange[ich][ieq][0] = rate_Trange[ich][jeq][0]
                    if rate_Trange[ich][ieq][0] == 0:
                        for jeq in range(len(rate_eq[ich])):
                            if rate_Trange[ich][jeq][0] != 0:
                                rate_Trange[ich][ieq][0] = rate_Trange[ich][jeq][0]
                for ieq in range(len(rate_eq[ich])-1):
                    if rate_Trange[ich][ieq][0] != rate_Trange[ich][ieq+1][0]:
                        rate_Trange[ich][ieq][0] = 4
                
                for ieq in range(len(rate_eq[ich])):
                    if 'for' in rate_eq[ich][ieq]:
                        Ts = re.findall('=(.*)~', rate_eq[ich][ieq])
                        Ts[0] = Ts[0].strip(' ')
                        if Ts[0] == '':
                            Ts[0] = '0'
                        rate_Trange[ich][ieq][1] = float(Ts[0])
                        Te = re.findall('~(.*)\[', rate_eq[ich][ieq])
                        Te[0] = Te[0].strip(' ')
                        if Te[0] == '':
                            Te[0] = '100000'
                        rate_Trange[ich][ieq][2] = float(Te[0])
                        if ieq >= 1:
                            if rate_Trange[ich][ieq-1][2] != rate_Trange[ich][ieq][1]:
                                action = 'error'
                                print('')
                                print('Please check the following reaction rate coefficient. Temperature ranges in each expression does not make sense.')
                                print('"for T = a ~ b [K]"')
                                print(' - Reaction: '+str(reaction_list[ich])+'')
                                print(' - Rate    : '+str(rate_eq[ich][ieq-1]).lstrip(' ')+' for T = '+str(rate_Trange[ich][ieq-1][1])+' ~ '+str(rate_Trange[ich][ieq-1][2])+' [K]')
                                print('           : '+str(rate_eq[ich][ieq]).lstrip(' ')+' for T = '+str(rate_Trange[ich][ieq][1])+' ~ '+str(rate_Trange[ich][ieq][2])+' [K]')

        ################# Converting rate expression into Reversed Polish Notation ################# 
        rate_cases[ich] = len(rate_eq[ich])
        for ieq in range(len(rate_eq[ich])): #number of T-dependence(coefficients) for each equation.

            #pressure dependent theree body reaction
            #it is recognized as pressure dependent three body reaction later in Fortran code
            if 'k0LH' in rate_eq[ich][ieq]:
                rate_Trange[ich][ieq][1] = -333
                rate_Trange[ich][ieq][2] = -333
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('k0LH','')
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('=','')
            if 'k2LH' in rate_eq[ich][ieq]:
                rate_Trange[ich][ieq][1] = -333
                rate_Trange[ich][ieq][2] = -333
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('k2LH','')
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('=','')
            if 'k3LH' in rate_eq[ich][ieq]:
                rate_Trange[ich][ieq][1] = -333
                rate_Trange[ich][ieq][2] = -333
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('k3LH','')
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('=','')
            if 'k0M' in rate_eq[ich][ieq]:
                rate_Trange[ich][ieq][1] = -333
                rate_Trange[ich][ieq][2] = -333
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('k0M','')
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('=','')
            if 'kinfM' in rate_eq[ich][ieq]: # coef4
                rate_Trange[ich][ieq][1] = -333
                rate_Trange[ich][ieq][2] = -333
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('kinfM','')
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('=','')
            if 'k0' in rate_eq[ich][ieq]:
                rate_Trange[ich][ieq][1] = -333
                rate_Trange[ich][ieq][2] = -333
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('k0','')
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('=','')
            if 'kinf' in rate_eq[ich][ieq]: # coef4
                rate_Trange[ich][ieq][1] = -333
                rate_Trange[ich][ieq][2] = -333
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('kinf','')
                rate_eq[ich][ieq] = rate_eq[ich][ieq].replace('=','')

            #special cases
            if 'Photoionization' in rate_eq[ich][ieq]:
                rate_eq[ich][ieq] = '0'
            if 'Photodissociation' in rate_eq[ich][ieq]:
                rate_eq[ich][ieq] = '0'
            if 'Impact' in rate_eq[ich][ieq]:
                rate_eq[ich][ieq] = '0'
            if 'Meteoroid' in rate_eq[ich][ieq]:
                rate_eq[ich][ieq] = '0'
            if 'Rainout' in rate_eq[ich][ieq]:
                rate_eq[ich][ieq] = '0'
            if 'for' in rate_eq[ich][ieq]:
                c = re.findall('(.*)for', rate_eq[ich][ieq])
                rate_eq[ich][ieq] = c[0]
            rate_eq[ich][ieq] = rate_eq[ich][ieq].lstrip()
            rate_eq[ich][ieq] = rate_eq[ich][ieq].rstrip()

            #### Convert Infix Notation to Reversed Polish Notation (RPN) ####
            Infix  = split_infix_tokens(rate_eq[ich][ieq])
            RPN    = infix_to_RPN(Infix)
            result = RPN_for_F90(RPN)

            rate_rpn_token[ich][ieq] = result[0]
            rate_rpn_label[ich][ieq] = result[1]

            if ntokens < rate_rpn_token[ich][ieq][0]:
                ntokens = rate_rpn_token[ich][ieq][0]

    ### For output rate information ###
    rate_rpn_token_for_f90  = [[[0 for ieq in range(ntokens+1)] for ieq in range(ncases)] for ich in range(len(rate_eq))]
    rate_rpn_label_for_f90  = [[[0 for ieq in range(ntokens+1)] for ieq in range(ncases)] for ich in range(len(rate_eq))]
    for ich in range(len(rate_rpn_token)):
        for ieq in range(len(rate_rpn_token[ich])):
            for itoken in range(len(rate_rpn_token[ich][ieq])):
                rate_rpn_token_for_f90[ich][ieq][itoken] = rate_rpn_token[ich][ieq][itoken]
                if rate_rpn_token_for_f90[ich][ieq][itoken] == 'T':
                    action = 'error'
                    print('')
                    print('Please indicate whether "T" is "Tn", "Ti" or "Te" in the following reaction rate coefficient.')
                    print('Tn: neutral temperature, Ti: ion temperature, Te: electron temperature.')
                    print(' - Reaction: '+str(reaction_list[ich])+'')
                    print(' - Rate    : '+str(rate_eq[ich][ieq])+'')
                rate_rpn_label_for_f90[ich][ieq][itoken] = rate_rpn_label[ich][ieq][itoken]

    # Output species info. on sub_window
    if action != 'error':
        species_info_win = tk.Toplevel()
        species_info_win.title("Species info")
        species_info_win.geometry("600x800")
        species_info = tk.Text(species_info_win, font=("",15), height=200, width=190, highlightthickness=0)
        species_info.pack()

        species_info.insert(tk.END, '############################################\n')
        species_info.insert(tk.END, '############### species info ###############\n')
        species_info.insert(tk.END, '############################################\n')
        species_info.insert(tk.END, '\n')
        species_info.insert(tk.END, '\n')
        nsp_f = 0
        for isp in range(len(species)):
            if fix_species[isp] == 1:
                nsp_f = nsp_f + 1
        species_info.insert(tk.END, '### Fixed : '+str(nsp_f)+' species ###\n')
        species_info.insert(tk.END, '\n')
        for isp in range(len(species)):
            species_for_info = species[isp]
            species_for_info = reaction_unicode(species_for_info)
            if fix_species[isp] == 1:
                species_info.insert(tk.END, '--------------------------------------------\n')
                species_info.insert(tk.END, ' '+species_for_info+'\n')
                species_info.insert(tk.END, '--------------------------------------------\n')
                for i in range(len(input_species_char)):
                    if species[isp] == input_species_char[i]:
                        path = input_species_path[i]
                        species_info.insert(tk.END, '# path    : '+path+'\n')
                species_info.insert(tk.END, '# mass    : '+'%.8f'%mass[isp]+'\n')
                species_info.insert(tk.END, '# charge  : '+charge[isp]+'\n')
                species_info.insert(tk.END, '\n')
        species_info.insert(tk.END, '--------------------------------------------\n')
        species_info.insert(tk.END, '\n')
        species_info.insert(tk.END, '\n')

        nsp_i = 0
        for isp in range(len(species)):
            if fix_species[isp] == 0:
                nsp_i = nsp_i + 1
        species_info.insert(tk.END, '### Variable : '+str(nsp_i)+' species ###\n')
        species_info.insert(tk.END, '\n')
        for isp in range(len(species)):
            species_for_info = species[isp]
            species_for_info = reaction_unicode(species_for_info)
            if fix_species[isp] == 0:
                species_info.insert(tk.END, '--------------------------------------------\n')
                species_info.insert(tk.END, ' '+species_for_info+'\n')
                species_info.insert(tk.END, '--------------------------------------------\n')

                for i in range(len(input_species_char)):
                    if species[isp] == input_species_char[i]:
                        path = input_species_path[i]
                        species_info.insert(tk.END, '# path    : '+path+'\n')
                species_info.insert(tk.END, '# mass    : '+'%.8f'%mass[isp]+'\n')
                species_info.insert(tk.END, '# charge  : '+charge[isp]+'\n')
                if Production_list[isp][0][0] != 0:
                    species_info.insert(tk.END, '# Produced by following '+str(Production_list[isp][0][0])+' reactions:\n')
                    species_info.insert(tk.END, '\n')
                    for ich in range(len(Production_list[isp])):
                        if ich != 0:
                            jch = Production_list[isp][ich][0]
                            if jch < 10:
                                species_info.insert(tk.END, '    * R'+str(jch)+'     '+reaction_list[jch-1]+'\n')
                            if jch >= 10 and jch < 100:
                                species_info.insert(tk.END, '    * R'+str(jch)+'    '+reaction_list[jch-1]+'\n')
                            if jch >= 100 and jch < 1000:
                                species_info.insert(tk.END, '    * R'+str(jch)+'   '+reaction_list[jch-1]+'\n')
                            if jch >= 1000 and jch < 10000:
                                species_info.insert(tk.END, '    * R'+str(jch)+'  '+reaction_list[jch-1]+'\n')
                    species_info.insert(tk.END, '\n')
                if Production_list[isp][0][0] == 0:
                    species_info.insert(tk.END, '# There are no reactions producing '+species_for_info+'.\n')
                    species_info.insert(tk.END, '\n')

                if Loss_list[isp][0][0] != 0:
                    species_info.insert(tk.END, '# Lost by following '+str(Loss_list[isp][0][0])+' reactions:\n')
                    species_info.insert(tk.END, '\n')
                    for ich in range(len(Loss_list[isp])):
                        if ich != 0:
                            jch = Loss_list[isp][ich][0]
                            if jch < 10:
                                species_info.insert(tk.END, '    * R'+str(jch)+'     '+reaction_list[jch-1]+'\n')
                            if jch >= 10 and jch < 100:
                                species_info.insert(tk.END, '    * R'+str(jch)+'    '+reaction_list[jch-1]+'\n')
                            if jch >= 100 and jch < 1000:
                                species_info.insert(tk.END, '    * R'+str(jch)+'   '+reaction_list[jch-1]+'\n')
                            if jch >= 1000 and jch < 10000:
                                species_info.insert(tk.END, '    * R'+str(jch)+'  '+reaction_list[jch-1]+'\n')
                    species_info.insert(tk.END, '\n')
                if Loss_list[isp][0][0] == 0:
                    species_info.insert(tk.END, '# There are no reactions losing '+species_for_info+'.\n')
                    species_info.insert(tk.END, '\n')
        species_info.insert(tk.END, '--------------------------------------------\n')

        # Output reaction info. on sub_window
        reaction_info_win = tk.Toplevel()
        reaction_info_win.title("reaction info")
        reaction_info_win.geometry("600x800")
        reaction_info = tk.Text(reaction_info_win, font=("",15), height=200, width=190, highlightthickness=0)
        reaction_info.pack()

        reaction_info.insert(tk.END, '############################################\n')
        reaction_info.insert(tk.END, '############## reaction info ###############\n')
        reaction_info.insert(tk.END, '############################################\n')
        reaction_info.insert(tk.END, '\n')
        reaction_info.insert(tk.END, '\n')
        reaction_info.insert(tk.END, '* '+str(len(reaction_list))+' reactions\n')
        reaction_info.insert(tk.END, '\n')

        for ich in range(len(reaction_list)):
            reaction_info.insert(tk.END, '--------------------------------------------\n')
            reaction_info.insert(tk.END, reaction_list[ich]+'\n')
            reaction_info.insert(tk.END, '--------------------------------------------\n')
            reaction_info.insert(tk.END, 'Reaction number : '+str(ich+1)+'\n')
            reaction_info.insert(tk.END, 'Reactant :  ')
            for isp in range(len(reactant[ich])):
                reactant_unicode = reaction_unicode(reactant[ich][isp])
                reaction_info.insert(tk.END, reactant_unicode)
                if isp != len(reactant[ich]) - 1:
                    reaction_info.insert(tk.END, '  ,  ')
            reaction_info.insert(tk.END, '\n')
            reaction_info.insert(tk.END, 'Product  :  ')
            for isp in range(len(product[ich])):
                product_unicode = reaction_unicode(product[ich][isp])
                reaction_info.insert(tk.END, product_unicode)
                if isp != len(product[ich]) - 1:
                    reaction_info.insert(tk.END, '  ,  ')
            reaction_info.insert(tk.END, '\n')
            reaction_info.insert(tk.END, '\n')
            if len(rate_eq[ich]) == 1:
                reaction_info.insert(tk.END, 'Reaction rate : '+str(len(rate_eq[ich]))+' case.\n')
            if len(rate_eq[ich]) != 1:
                reaction_info.insert(tk.END, 'Reaction rate : '+str(len(rate_eq[ich]))+' cases.\n')
            for ieq in range(len(rate_eq[ich])):
                rate_eq_unicode = rate_unicode(rate_eq[ich][ieq])
                if rate_Trange[ich][ieq][1] >= 0:
                    reaction_info.insert(tk.END, '       ('+str(ieq+1)+') '+rate_eq_unicode +'\n')
                if rate_Trange[ich][ieq][1] == -333 and rate_Trange[ich][ieq][2] == -333:
                    reaction_info.insert(tk.END, '       ('+str(ieq+1)+') k0 = '+rate_eq_unicode +'\n')
                if rate_Trange[ich][ieq][1] == -333 and rate_Trange[ich][ieq][2] > 1000:
                    reaction_info.insert(tk.END, '       ('+str(ieq+1)+') kinf = '+rate_eq_unicode +'\n')
                #reaction_info.insert(tk.END, '               type  :  ')
                #if rate_label[ich][ieq] == 0:
                #    reaction_info.insert(tk.END, 'a (constant)\n')
                #    if rate_coef[ich][ieq][0] > 0:
                #        reaction_info.insert(tk.END, '                 a = '+'%9.5e'%rate_coef[ich][ieq][1]+'\n')
                #    if rate_coef[ich][ieq][0] == -1:
                #        reaction_info.insert(tk.END, '                 Photoionization\n')
                #    if rate_coef[ich][ieq][0] == -2:
                #        reaction_info.insert(tk.END, '                 Impact ionization\n')
                #    if rate_coef[ich][ieq][0] == -3:
                #        reaction_info.insert(tk.END, '                 Meteoroid ablation\n')
                #    if rate_coef[ich][ieq][0] == -4:
                #        reaction_info.insert(tk.END, '                 Condensation\n')
                #if rate_coef[ich][ieq][0] == 1:
                #    reaction_info.insert(tk.END, 'a \u00D7 (b/T)^c\n')
                #    reaction_info.insert(tk.END, '                 a = '+'%6.2e'%rate_coef[ich][ieq][1]+'\n')
                #    reaction_info.insert(tk.END, '                 b = '+'%6.2f'%rate_coef[ich][ieq][2]+'\n')
                #    reaction_info.insert(tk.END, '                 c = '+'%3.2f'%rate_coef[ich][ieq][3]+'\n')
                #if rate_coef[ich][ieq][0] == 2:
                #    reaction_info.insert(tk.END, 'a \u00D7 exp(b/T)\n')
                #    reaction_info.insert(tk.END, '                 a = '+'%6.2e'%rate_coef[ich][ieq][1]+'\n')
                #    reaction_info.insert(tk.END, '                 b = '+'%6.2f'%rate_coef[ich][ieq][2]+'\n')
                #if rate_coef[ich][ieq][0] == 3:
                #    reaction_info.insert(tk.END, 'a \u00D7 (b/T)^c \u00D7 exp(d/T)\n')
                #    reaction_info.insert(tk.END, '                 a = '+'%6.2e'%rate_coef[ich][ieq][1]+'\n')
                #    reaction_info.insert(tk.END, '                 b = '+'%6.2f'%rate_coef[ich][ieq][2]+'\n')
                #    reaction_info.insert(tk.END, '                 c = '+'%3.2f'%rate_coef[ich][ieq][3]+'\n')
                #    reaction_info.insert(tk.END, '                 d = '+'%6.2f'%rate_coef[ich][ieq][4]+'\n')
                #if rate_coef[ich][ieq][0] == 4:
                #    reaction_info.insert(tk.END, 'Three body reaction\n')
                #    #reaction_info.insert(tk.END, 'a \u00D7 (b/T)^c \u00D7 exp(d/T)\n')
                #    #reaction_info.insert(tk.END, '                 a = '+'%6.2e'%rate_coef3[ich][ieq][0]+'\n')
                #    #reaction_info.insert(tk.END, '                 b = '+'%6.2f'%rate_coef3[ich][ieq][1]+'\n')
                #    #reaction_info.insert(tk.END, '                 c = '+'%3.2f'%rate_coef3[ich][ieq][2]+'\n')
                #    #reaction_info.insert(tk.END, '                 d = '+'%6.2f'%rate_coef3[ich][ieq][3]+'\n')
                reaction_info.insert(tk.END, '               T type  : ')
                if rate_Trange[ich][ieq][0] == 0:
                    reaction_info.insert(tk.END, '\n')
                if rate_Trange[ich][ieq][0] == 1:
                    reaction_info.insert(tk.END, 'neutral Temperature\n')
                if rate_Trange[ich][ieq][0] == 2:
                    reaction_info.insert(tk.END, 'electron Temperature\n')
                if rate_Trange[ich][ieq][0] == 3:
                    reaction_info.insert(tk.END, 'ion Temperature\n')
                if rate_Trange[ich][ieq][0] == 4:
                    reaction_info.insert(tk.END, 'error\n')
                    print('error : '+reaction_list[ich]+'  '+rate_eq[ich][ieq])
                reaction_info.insert(tk.END, '               T range : ')
                if rate_Trange[ich][ieq][1] >= 0:
                    reaction_info.insert(tk.END, '%6.2f'%rate_Trange[ich][ieq][1]+' [K]  ~  '+'%6.2f'%rate_Trange[ich][ieq][2]+' [K]\n')
                if rate_Trange[ich][ieq][1] < 0:
                    reaction_info.insert(tk.END, '0.00 [K]  ~  100000.00 [K]\n')
                reaction_info.insert(tk.END, '\n')

            reaction_info.insert(tk.END, '\n')
            reaction_info.insert(tk.END, '\n')


    ##############################################################################################################
    ##################################   Output Photochemical module   ###########################################
    ##############################################################################################################
    if action == 'Output' or action == 'Run':

        #### create necessary directories ########################################
        path = './build'
        if os.path.exists(path) == False:
            os.makedirs(path)

        path = './UV/xct/absorption'
        if os.path.exists(path) == False:
            os.makedirs(path)
        
        path = './UV/xct/dissociation'
        if os.path.exists(path) == False:
            os.makedirs(path)

        path = './UV/xct/RS'
        if os.path.exists(path) == False:
            os.makedirs(path)

        #### output ################################################################
        #write selected reaction list
        path = './'+Planet+'/'+dir0+'/settings/reaction_list.dat'
        with open(path, mode = 'w') as f:
            for i in range(len(reaction_rate_list)):
                if reaction_chk_bln[i].get() == True: #only checked reaction
                    f.write(reaction_rate_list[i]+'\n')

        path = './'+Planet+'/'+dir0+'/settings/reaction_list_for_paper.csv'
        with open(path, mode = 'w') as f:
            for i in range(len(reaction_list)):
                char = reaction_unicode(reaction_list_L[i])
                f.write(char+', \u2192 ,')
                char = reaction_unicode(reaction_list_R[i])
                f.write(char+',')
                for j in range(len(rate_eq_list[i])):
                    char = re.sub('\s+', ' ', rate_eq_list[i][j])
                    char = rate_unicode(char)
                    if j == 0:
                        f.write(char+','+ref_label[i][j]+'\n')
                    if j >= 1:
                        f.write(' , , ,'+char+','+ref_label[i][j]+'\n')

        path = './'+Planet+'/'+dir0+'/settings/reaction_list_for_paper.tex'
        with open(path, mode = 'w') as f:
          
            f.write('\\documentclass[10pt,a4paper]{jarticle}\n')
            f.write('\\usepackage[dvipdfmx]{graphicx}\n')
            f.write('\\usepackage{geometry}\n')
            f.write('\\usepackage{here}\n')
            f.write('\\usepackage{graphicx}\n')
            f.write('\\usepackage{tabularx}\n')
            f.write('\\usepackage{scalefnt}\n')
            f.write('\\geometry{left=10mm,right=25mm,top=30mm,bottom=30mm}\n')
            f.write('\n')
            f.write('\\begin{document}\n')
            f.write('\n')
            f.write('\\begin{table}[htb]\n')
            f.write('\\scalebox{0.7}{')
            f.write('\\begin{tabular}{lclll} \\hline\n')
            f.write('    \\multicolumn{3}{c}{Reactions} & \\multicolumn{1}{c}{Rate [cm$^{3}$/s]} & \\multicolumn{1}{c}{Reference}\\\\\\hline\n')
            
            for i in range(len(reaction_list)):
                char = reaction_unicode(reaction_list_L[i])
                char = unicode_LaTeX(char)
                f.write(char+'& $\\rightarrow$ &')
                char = reaction_unicode(reaction_list_R[i])
                char = unicode_LaTeX(char)
                f.write(char+'&')
                for j in range(len(rate_eq_list[i])):
                    char = re.sub('\s+', ' ', rate_eq_list[i][j])
                    char = rate_unicode(char)
                    char = unicode_LaTeX(char)
                    if j == 0:
                        f.write(char+'&'+ref_label[i][j]+'\\\\\n')
                    if j >= 1:
                        f.write(' & & &'+char+'&'+ref_label[i][j]+'\\\\\n')
            f.write('    \\end{tabular}\n')
            f.write('    }\n')
            f.write('\\end{table}\n')
            f.write('\\end{document}\n')


        #write reaction type as text file
        path = './'+Planet+'/'+dir0+'/input/PLJ_list/reaction_type_list.dat'
        with open(path, mode = 'w') as f:
            for ich in range(len(reaction_type_list)):
                f.write(str(reaction_type_list[ich])+' ')
                f.write('\n')

        path = './'+Planet+'/'+dir0+'/input/PLJ_list/reactant_list.dat'
        with open(path, mode = 'w') as f:
            for ich in range(len(reactant_list)):
                for isp in range(len(reactant_list[ich])):
                    f.write(str(reactant_list[ich][isp])+' ')
                f.write('\n')

        path = './'+Planet+'/'+dir0+'/input/PLJ_list/product_list.dat'
        with open(path, mode = 'w') as f:
            for ich in range(len(product_list)):
                for isp in range(len(product_list[ich])):
                    f.write(str(product_list[ich][isp])+' ')
                f.write('\n')

        # write Production list
        path = './'+Planet+'/'+dir0+'/input/PLJ_list/Production_list.dat'
        with open(path, mode = 'w') as f:
            for isp in range(len(species)):
                if fix_species[isp] == 0:
                    for ich in range(len(Production_list_for_output[isp])):
                        f.write(str(Production_list_for_output[isp][ich])+' ')
                    f.write('\n')

        # write Loss list
        path = './'+Planet+'/'+dir0+'/input/PLJ_list/Loss_list.dat'
        with open(path, mode = 'w') as f:
            for isp in range(len(species)):
                if fix_species[isp] == 0:
                    for ich in range(len(Loss_list_for_output[isp])):
                        f.write(str(Loss_list_for_output[isp][ich])+' ')
                    f.write('\n')

        # write Jacbobian list
        path = './'+Planet+'/'+dir0+'/input/PLJ_list/Jacobian_list.dat'
        with open(path, mode = 'w') as f:
            # Jacobian_label is a matrix composed of 1 or 0
            # 0-> non-variable
            # 1-> has variable
            # only Jacobian_label(i,j) of 1 is read in Fortran.
            # this is needed to run faster!

            #the number of (i,j) of 1
            n_Jlist = 0

            # make new Jacobian only composed of variable species(without fixed species)
            # i0 and j0 is new index of Jacobian without fixed_species
            i0 = 0
            for isp in range(len(species)):
                if fix_species[isp] == 0:
                    i0 = i0 + 1
                    j0 = 0
                    for jsp in range(len(species)):
                        if fix_species[jsp] == 0:
                            j0 = j0 + 1
                            if Jacobian_label_P[isp][jsp] == 1 or Jacobian_label_L[isp][jsp] == 1: #if this (i0,j0) has some contribution(if not 0)
                                n_Jlist = n_Jlist + 1
                                f.write(str(i0)+' '+str(j0)+' ')
                                for i in range(max_ch_J):
                                    f.write(str(Jacobian_list_for_output[isp][jsp][i])+' ')
                                f.write('\n')


        path = './'+Planet+'/'+dir0+'/input/PLJ_list/rate_rpn_token.dat'
        with open(path, mode = 'w') as f:
            for ich in range(len(rate_rpn_token_for_f90)):
                for ieq in range(len(rate_rpn_token_for_f90[ich])):
                    for itoken in range(len(rate_rpn_token_for_f90[ich][ieq])):
                        f.write(str(rate_rpn_token_for_f90[ich][ieq][itoken])+' ')
                    f.write('\n')

        path = './'+Planet+'/'+dir0+'/input/PLJ_list/rate_rpn_label.dat'
        with open(path, mode = 'w') as f:
            for ich in range(len(rate_rpn_label_for_f90)):
                for ieq in range(len(rate_rpn_label_for_f90[ich])):
                    for itoken in range(len(rate_rpn_label_for_f90[ich][ieq])):
                        f.write(str(rate_rpn_label_for_f90[ich][ieq][itoken])+' ')
                    f.write('\n')

        path = './'+Planet+'/'+dir0+'/input/PLJ_list/rate_cases.dat'
        with open(path, mode = 'w') as f:
            for ich in range(len(rate_eq)):
                f.write(str(len(rate_eq[ich]))+' ')
            f.write('\n')

        path = './'+Planet+'/'+dir0+'/input/PLJ_list/T_range.dat'
        with open(path, mode = 'w') as f:
            for ich in range(len(rate_Trange)):
                for ieq in range(len(rate_Trange[ich])):
                    for i in range(len(rate_Trange[ich][ieq])):
                        f.write(str(rate_Trange[ich][ieq][i])+' ')
                    f.write('\n')


        # Output species info.
        path = './'+Planet+'/'+dir0+'/info/species_info.dat'
        with open(path, mode = 'w') as f:
            f.write(species_info.get('1.0', tk.END))
        f.close()
        # Output reaction info.
        path = './'+Planet+'/'+dir0+'/info/reaction_info.dat'
        with open(path, mode = 'w') as f:
            f.write(reaction_info.get('1.0', tk.END))
        f.close()

        # read grid data
        dz = {}
        zrange = {}
        path = './'+Planet+'/'+dir0+'/settings/zgrid.dat'
        if os.path.exists(path) == True:
            with open(path, mode = 'r') as f:
                lines = f.readlines()
                for iz in range(len(lines)):
                    lchar = lines[iz]
                    rchar = lines[iz]
                    lchar = lchar.split(',')
                    dz[iz] = lchar[0]
                    rchar = rchar.split(',')
                    rchar = rchar[1].split('\n')
                    zrange[iz] = rchar[0].lstrip()
            f.close()
            if len(lines)==0:
                err_win = tk.Toplevel() #display on the main window
                err_win.title("error")
                err_win.geometry("300x100")
                txt='Please set vertical grid!'
                txt_label = tk.Label(err_win, anchor="w", justify="left", font=('', 15), text=txt)
                txt_label.place(x=10, y=10)
        if os.path.exists(path) == False:
            err_win = tk.Toplevel() #display on the main window
            err_win.title("error")
            err_win.geometry("300x100")
            txt='Please set vertical grid!'
            txt_label = tk.Label(err_win, anchor="w", justify="left", font=('', 15), text=txt)
            txt_label.place(x=10, y=10)

        nalt = 0
        nalt_iz = {}
        zs = {}
        ze = {}
        for iz in range(len(dz)):
            char = zrange[iz].split(':')
            zs[iz] = float(char[0])
            ze[iz] = float(char[1])
            nalt_iz[iz,1] = nalt + 1
            nalt_iz[iz,2] = nalt + (ze[iz]-zs[iz])/float(dz[iz])
            nalt = nalt + (ze[iz]-zs[iz])/float(dz[iz]) + 1

        path = './'+Planet+'/'+dir0+'/settings/calculation_setting.dat'
        dimension_input = ''
        lat_input = ''
        f107_input = ''
        ls_input = ''
        stepmax_input = ''
        sumdt_input = ''
        dtmax_input = ''
        sumdtrot_input = ''
        dtmaxrot_input = ''
        scheme_input = ''
        inversion_input = ''
        if os.path.exists(path) == True:
            with open(path, mode = 'r') as f:
                line = f.readlines()
            f.close()

            for i in range(len(line)):
                if 'Dimension' in line[i]:
                    dimension_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
                if 'LAT' in line[i]:
                    lat_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
                if 'F10.7' in line[i]:
                    f107_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
                if 'Ls' in line[i]:
                    ls_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
                if 'stepmax' in line[i]:
                    stepmax_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
                if 'sumdtst' in line[i]:
                    sumdt_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
                if 'dtmaxst' in line[i]:
                    dtmax_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
                if 'sumdtrot' in line[i]:
                    sumdtrot_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
                if 'dtmaxrot' in line[i]:
                    dtmaxrot_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
                if 'scheme' in line[i]:
                    scheme_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
                if 'inversion' in line[i]:
                    inversion_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()


        ##### save as 'v__in.f90'
        print('Output f90 module...')
        path = './'+Planet+'/'+dir0+'/v__in.f90'
        with open(path, mode = 'w') as f:

            f.write(   'module v__in\n')
            f.write(   '  use v__tdec,   only : set_, var_, cst_, spl_, grd_\n')
            f.write(   '  use p__search, only : p__search_reactant, p__search_product, sp_index\n')
            f.write(   '\n')
            f.write(   '  implicit none\n')
            f.write(   '  integer(4), parameter         :: sp = 4, dp = 8\n')
            f.write(   '\n')
            f.write(   '  private\n')
            f.write(   '  public :: v__in__ini, v__in__exe\n')
            f.write(   '\n')
            f.write(   'contains\n')
            f.write(   '\n')
            f.write(   '\n')
            f.write(   '  subroutine v__in__ini(spl, set) ! out\n')
            f.write(   '    type(spl_),           intent(out) :: spl\n')
            f.write(   '    type(set_),           intent(out) :: set\n')
            f.write(   '\n')
            f.write(   '    ! Planet type\n')
            f.write(   '    spl%planet = \''+Planet+'\'\n')
            f.write(   '\n')
            f.write(   '    ! Calculation settings\n')
            dimension_char = ''
            if dimension_input == '1':
                dimension_char = '1D'
            if dimension_input == '2':
                dimension_char = '2D Lat'
            if dimension_input == '3':
                dimension_char = '2D Rot'
            if dimension_input == '4':
                dimension_char = '3D Rot'
            f.write(   '    set%mode = \''+dimension_char+'\'\n')
            f.write(   '    set%F107 = '+f107_input+'_dp\n')
            f.write(   '    set%nstep = '+stepmax_input+'\n')
            f.write(   '    set%fin_sec = '+sumdt_input+'_dp\n')
            f.write(   '    set%dtime_limit = '+dtmax_input+'_dp\n')
            f.write(   '    set%latitude = '+lat_input+'_dp\n')
            f.write(   '    set%Ls = '+ls_input+'_dp\n')
            f.write(   '    set%nday = '+sumdtrot_input+'_dp\n')
            scheme_char = ''
            if scheme_input == '1':
                scheme_char = 'implicit'
            if scheme_input == '2':
                scheme_char = 'semi-implicit'
            if scheme_input == '3':
                scheme_char = 'explicit'
            f.write(   '    set%scheme = \''+scheme_char+'\'\n')
            inv_char = ''
            if inversion_input == '1':
                inv_char = 'default'
            f.write(   '    set%inversion = \''+inv_char+'\'\n')
            f.write(   '    ! directory setting\n')
            f.write(   '    set%dir_name = \'./'+Planet+'/'+dir0+'\'\n')
            f.write(   '\n')
            f.write(   '  end subroutine v__in__ini\n')
            f.write(   '\n')
            f.write(   '\n')
            f.write(   '  subroutine v__in__exe(cst,      & ! in\n')
            f.write(   '    &                   set, spl, & ! inout\n')
            f.write(   '    &                   var, grd  ) ! out\n')
            f.write(   '    type(cst_),           intent(in)    :: cst\n')
            f.write(   '    type(set_),           intent(inout) :: set\n')
            f.write(   '    type(spl_),           intent(inout) :: spl\n')
            f.write(   '    type(var_),           intent(out)   :: var\n')
            f.write(   '    type(grd_),           intent(out)   :: grd\n')
            f.write(   '    integer i, j, ip, isp, jsp, ich, jch, iz, nh\n')
            f.write(   '    real(dp) tmp\n')
            f.write(   '    character(len = 256) strm, fname\n')
            f.write(   '\n')
            f.write(   '\n')

            #f.write(   '    integer, parameter :: zs = '+str(zs)+'\n')
            f.write(   '    ! grid setting\n')
            f.write(   '    grd%nx    = set%nx\n')
            f.write(   '    grd%ny    = set%ny\n')
            f.write(   '    grd%nz    = '+str(int(nalt))+'\n')
            f.write(   '    allocate(grd%alt(grd%nz),grd%dalt(grd%nz))\n')
            for iz in range(len(dz)):
                f.write(   '    grd%dalt('+str(int(nalt_iz[iz,1]))+':'+str(int(nalt_iz[iz,2])+1)+') = '+str(dz[iz])+'e3_dp ! [m]\n')
            f.write(   '    grd%alt(1)      = '+str(zs[0])+'e3_dp ! [m]\n')
            f.write(   '    do iz = 2, grd%nz\n')
            f.write(   '      grd%alt(iz) = grd%alt(iz-1) + grd%dalt(iz)\n')
            f.write(   '    end do\n')
            f.write(   '\n')
            f.write(   '    ! reactions, chemical species\n')
            f.write(   '    spl%nsp     = '+str(nsp_i+nsp_f)+'\n')
            f.write(   '    spl%nsp_i   = '+str(nsp_i)+'\n')
            f.write(   '    spl%nch     = '+str(len(reactant))+'\n')
            f.write(   '    spl%nch_P   = '+str(max_ch_P)+'\n')
            f.write(   '    spl%nch_L   = '+str(max_ch_L)+'\n')
            f.write(   '    spl%n_Jlist = '+str(n_Jlist)+'\n')
            f.write(   '    spl%nch_J   = '+str(max_ch_J+2)+'\n')
            f.write(   '    spl%nrpn    = '+str(ntokens+1)+'\n')
            f.write(   '\n')
            f.write(   '    ! allocate\n')
            f.write(   '    allocate(var%ni(0:spl%nsp,grd%nz))\n')
            f.write(   '    allocate(var%ni_0(0:spl%nsp,grd%nz))\n')
            f.write(   '    allocate(var%ni_new(spl%nsp,grd%nz))\n')
            f.write(   '    allocate(var%ni_stable(spl%nsp,grd%ny,grd%nz))\n')
            f.write(   '    allocate(var%ni_3d(spl%nsp,grd%nx,grd%ny,grd%nz))\n')
            f.write(   '    allocate(var%clm_ni(spl%nsp,grd%nz))\n')
            f.write(   '    allocate(var%Ti(grd%nz),var%Te(grd%nz),var%Tn(grd%nz))\n')
            f.write(   '    allocate(var%Ti_3d(grd%nx,grd%ny,grd%nz))\n')
            f.write(   '    allocate(var%Te_3d(grd%nx,grd%ny,grd%nz))\n')
            f.write(   '    allocate(var%Tn_3d(grd%nx,grd%ny,grd%nz))\n')
            f.write(   '    allocate(var%m(spl%nsp), var%q(spl%nsp))\n')
            f.write(   '    allocate(spl%reactant_list(spl%nch, 10))\n')
            f.write(   '    allocate(spl%product_list(spl%nch, 10))\n')
            f.write(   '    allocate(spl%species(0:spl%nsp))\n')
            f.write(   '    allocate(spl%label_fix(spl%nsp))\n')
            f.write(   '    allocate(spl%all_to_var(spl%nsp))\n')
            f.write(   '    allocate(spl%var_to_all(spl%nsp_i))\n')
            f.write(   '    allocate(spl%Prod_list(spl%nsp_i, spl%nch_P))\n')
            f.write(   '    allocate(spl%Loss_list(spl%nsp_i, spl%nch_L))\n')
            f.write(   '    allocate(spl%Jmtx_list(spl%n_Jlist, spl%nch_J))\n')
            f.write(   '    allocate(spl%reaction_type_list(spl%nch))\n')
            f.write(   '    allocate(spl%reaction_type_char(spl%nch))\n')
            f.write(   '    allocate(var%ki(spl%nch,grd%nz))\n')
            f.write(   '    allocate(var%rate(spl%nch,grd%nz))\n')
            f.write(   '    allocate(var%Pi(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%Pij(spl%nsp_i,grd%nz,spl%nch))\n')
            f.write(   '    allocate(var%Li(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%K_eddy(grd%nz),var%D_mol(spl%nsp,grd%nz))\n')
            f.write(   '    allocate(var%Fluxup(spl%nsp_i,0:grd%nz),var%Fluxdwn(spl%nsp_i,0:grd%nz))\n')
            f.write(   '    allocate(var%vFluxup(spl%nsp_i,grd%nz),var%vFluxdwn(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%UpperBC(0:spl%nsp,3),var%LowerBC(0:spl%nsp,3))\n')
            f.write(   '    allocate(var%Jmtx(spl%nsp_i, spl%nsp_i))\n')
            f.write(   '    allocate(spl%rate_rpn_token(spl%nch,3,spl%nrpn))\n')
            f.write(   '    allocate(spl%rate_rpn_label(spl%nch,3,spl%nrpn))\n')
            f.write(   '    allocate(spl%rate_cases(spl%nch))\n')
            f.write(   '    allocate(spl%T_range(spl%nch,3,3))\n')
            f.write(   '    allocate(spl%major_species(grd%nz))\n')
            f.write(   '    allocate(var%n_tot(grd%nz),var%m_mean(grd%nz))\n')
            f.write(   '    allocate(var%Phip(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%Phim(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%dPhi_dz(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%d_dniu_dPhi_dz(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%d_dni0_dPhi_dz(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%d_dnil_dPhi_dz(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%d_dneu_dPhi_dz_add(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%d_dne0_dPhi_dz_add(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%d_dnel_dPhi_dz_add(spl%nsp_i,grd%nz))\n')
            f.write(   '    allocate(var%barr(spl%nsp_i*grd%nz), var%xarr(spl%nsp_i*grd%nz))\n')
            f.write(   '    allocate(var%yarr(spl%nsp_i*grd%nz), var%dxarr(spl%nsp_i*grd%nz))\n')
            f.write(   '    allocate(var%tAmtx(2*spl%nsp_i+1,spl%nsp_i*grd%nz))\n')
            f.write(   '    allocate(var%tLmtx(spl%nsp_i+1,spl%nsp_i*grd%nz))\n')
            f.write(   '    allocate(var%Umtx(spl%nsp_i+1,spl%nsp_i*grd%nz))\n')

            f.write(   '\n')
            f.write(   '    ! species\n')
            f.write(   '    ! \n')
            f.write(   '    ! ')
            for isp in range(len(species)-1):
                f.write(species[isp]+', ')
            f.write(species[len(species)-1]+'\n')
            f.write(   '    ! \n')
            jsp = 0
            for isp in range(len(species)):
                #if fix_species[isp] == 0:
                f.write(   '    spl%species('+str(jsp+1)+') = \''+species[isp]+'\'\n')
                jsp = jsp + 1

            f.write(   '\n')
            f.write(   '    ! label_fix\n')
            jsp = 0
            for isp in range(len(species)):
                if fix_species[isp] == 0:
                    f.write(   '    spl%label_fix('+str(jsp+1)+') = 0 ! '+species[isp]+': variable\n')
                    jsp = jsp + 1
                if fix_species[isp] == 1:
                    f.write(   '    spl%label_fix('+str(jsp+1)+') = 1 ! '+species[isp]+': fixed\n')
                    jsp = jsp + 1

            f.write(   '\n')
            f.write(   '    ! all_to_var\n')
            f.write(   '    spl%all_to_var = 0\n')
            jsp = 0
            for isp in range(len(species)):
                if fix_species[isp] == 0:
                    f.write(   '    spl%all_to_var('+str(isp+1)+') = '+str(jsp+1)+' ! '+species[isp]+': variable\n')
                    jsp = jsp + 1
            f.write(   '\n')
            f.write(   '    ! var_to_all\n')
            jsp = 0
            for isp in range(len(species)):
                if fix_species[isp] == 0:
                    f.write(   '    spl%var_to_all('+str(jsp+1)+') = '+str(isp+1)+' ! '+species[isp]+': variable\n')
                    jsp = jsp + 1

            f.write(   '\n')
            f.write(   '    ! mass\n')
            jsp = 0
            for isp in range(len(species)):
                #if fix_species[isp] == 0:
                f.write(   '    var%m('+str(jsp+1)+') = '+'%.8f'%mass[isp]+'_dp * cst%m_u !'+species[isp]+'\n')
                jsp = jsp + 1
            f.write(   '\n')
            f.write(   '    ! mass zero error\n')
            f.write(   '    do isp = 1, spl%nsp\n')
            f.write(   '      if ( var%m(isp) == 0.0_dp ) then \n')
            f.write(   '        write(*,*) \'mass zero error!\'\n')
            f.write(   '        write(*,*) \'mass of \',trim(ADJUSTL(spl%species(isp))),\' is zero.\'\n')
            f.write(   '        write(*,*) \'Calculation stopped.\'\n')
            f.write(   '        stop\n')
            f.write(   '      end  if\n')
            f.write(   '    end do\n')

            f.write(   '\n')
            f.write(   '    ! charge\n')
            jsp = 0
            for isp in range(len(species)):
                #if fix_species[isp] == 0:
                f.write(   '    var%q('+str(jsp+1)+') = '+str(charge[isp])+'.0_dp * cst%q_e !'+species[isp]+'\n')
                jsp = jsp + 1

            f.write(   '\n')
            f.write(   '    ! read P, L, J list\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/PLJ_list/Production_list.dat\', status = \'unknown\' )\n')
            f.write(   '      do isp = 1, spl%nsp_i\n')
            f.write(   '        read(11,*) (spl%Prod_list(isp,ich), ich = 1, spl%nch_P)\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/PLJ_list/Loss_list.dat\', status = \'unknown\' )\n')
            f.write(   '      do isp = 1, spl%nsp_i\n')
            f.write(   '        read(11,*) (spl%Loss_list(isp,ich), ich = 1, spl%nch_L)\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/PLJ_list/Jacobian_list.dat\', status = \'unknown\' )\n')
            f.write(   '      do i = 1, spl%n_Jlist\n')
            f.write(   '        read(11,*) (spl%Jmtx_list(i,j), j = 1, spl%nch_J)\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/PLJ_list/reactant_list.dat\', status = \'unknown\' )\n')
            f.write(   '      do ich = 1, spl%nch\n')
            f.write(   '        read(11,*) (spl%reactant_list(ich,isp), isp = 1, 10)\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/PLJ_list/product_list.dat\', status = \'unknown\' )\n')
            f.write(   '      do ich = 1, spl%nch\n')
            f.write(   '        read(11,*) (spl%product_list(ich,isp), isp = 1, 10)\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/PLJ_list/reaction_type_list.dat\', status = \'unknown\' )\n')
            f.write(   '      do ich = 1, spl%nch\n')
            f.write(   '        read(11,*) spl%reaction_type_list(ich)\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/PLJ_list/rate_rpn_token.dat\', status = \'unknown\' )\n')
            f.write(   '      do ich = 1, spl%nch\n')
            f.write(   '        do i = 1, 3\n')
            f.write(   '          read(11,*) (spl%rate_rpn_token(ich,i,j), j = 1, spl%nrpn)\n')
            f.write(   '        end do\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/PLJ_list/rate_rpn_label.dat\', status = \'unknown\' )\n')
            f.write(   '      do ich = 1, spl%nch\n')
            f.write(   '        do i = 1, 3\n')
            f.write(   '          read(11,*) (spl%rate_rpn_label(ich,i,j), j = 1, spl%nrpn)\n')
            f.write(   '        end do\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/PLJ_list/rate_cases.dat\', status = \'unknown\' )\n')
            f.write(   '      read(11,*) (spl%rate_cases(ich), ich = 1, spl%nch)\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/PLJ_list/T_range.dat\', status = \'unknown\' )\n')
            f.write(   '      do ich = 1, spl%nch\n')
            f.write(   '        do i = 1, 3\n')
            f.write(   '          read(11,*) (spl%T_range(ich,i,j), j = 1, 3)\n')
            f.write(   '        end do\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    ! reaction type characters\n')
            f.write(   '    var%nspecial = 0\n')
            f.write(   '    do ich = 1, spl%nch\n')
            f.write(   '      if (      spl%reaction_type_list(ich) == '+str(type_photoionization)+' ) then\n')
            f.write(   '        spl%reaction_type_char(ich) = \'photoionization\'\n')
            f.write(   '      else if ( spl%reaction_type_list(ich) == '+str(type_photodissociaion)+' ) then\n')
            f.write(   '        spl%reaction_type_char(ich) = \'photodissociation\'\n')
            f.write(   '      else if ( spl%reaction_type_list(ich) == '+str(type_e_impact)+' ) then\n')
            f.write(   '        spl%reaction_type_char(ich) = \'electron impact\'\n')
            f.write(   '        var%nspecial = var%nspecial + 1\n')
            f.write(   '      else if ( spl%reaction_type_list(ich) == '+str(type_p_impact)+' ) then\n')
            f.write(   '        spl%reaction_type_char(ich) = \'proton impact\'\n')
            f.write(   '        var%nspecial = var%nspecial + 1\n')
            f.write(   '      else if ( spl%reaction_type_list(ich) == '+str(type_H_impact)+' ) then\n')
            f.write(   '        spl%reaction_type_char(ich) = \'H impact\'\n')
            f.write(   '        var%nspecial = var%nspecial + 1\n')
            f.write(   '      else if ( spl%reaction_type_list(ich) == '+str(type_pressure_3body)+' ) then\n')
            f.write(   '        spl%reaction_type_char(ich) = \'pressure_dependent_3body\'\n')
            f.write(   '      else if ( spl%reaction_type_list(ich) == '+str(type_pressure_3bodyM)+' ) then\n')
            f.write(   '        spl%reaction_type_char(ich) = \'pressure_dependent_3bodyM\'\n')
            f.write(   '      else if ( spl%reaction_type_list(ich) == '+str(type_Lindemann_Hinshelwood)+' ) then\n')
            f.write(   '        spl%reaction_type_char(ich) = \'Lindemann-Hinshelwood\'\n')
            f.write(   '      else if ( spl%reaction_type_list(ich) == '+str(type_Meteoroid)+' ) then\n')
            f.write(   '        spl%reaction_type_char(ich) = \'Meteoroid ablation\'\n')
            f.write(   '        var%nspecial = var%nspecial + 1\n')
            f.write(   '      else if ( spl%reaction_type_list(ich) == '+str(type_Rainout)+' ) then\n')
            f.write(   '        spl%reaction_type_char(ich) = \'Rainout\'\n')
            f.write(   '        var%nspecial = var%nspecial + 1\n')
            f.write(   '      end if\n')
            f.write(   '    end do\n')

            f.write(   '\n')

            # temperature input
            f.write(   '    ! input Temperature profiles\n')
            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/Temperature/T_e.dat\', status = \'unknown\' )\n')
            f.write(   '      do iz = 1, grd%nz\n')
            f.write(   '        read(11,*) tmp,var%Te(iz)\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/Temperature/T_i.dat\', status = \'unknown\' )\n')
            f.write(   '      do iz = 1, grd%nz\n')
            f.write(   '        read(11,*) tmp,var%Ti(iz)\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            f.write(   '\n')
            f.write(   '    open(11, file = \'./'+Planet+'/'+dir0+'/input/Temperature/T_n.dat\', status = \'unknown\' )\n')
            f.write(   '      do iz = 1, grd%nz\n')
            f.write(   '        read(11,*) tmp,var%Tn(iz)\n')
            f.write(   '      end do\n')
            f.write(   '    close(11)\n')

            # neutral input
            f.write(   '\n')
            f.write(   '    ! input density profiles\n')
            f.write(   '    var%ni   = 1.0e-20_dp\n')
            for isp in range(len(species)):
                for i in range(len(input_species_char)):
                    if species[isp] == input_species_char[i]:
                        f.write(   '\n')
                        f.write(   "    isp = sp_index(spl, '"+species[isp]+"')\n")
                        f.write(   '    if (isp >= 1 .and. isp <= spl%nsp) then\n')
                        f.write(   '      open(11, file = \''+input_species_path[i]+'\', status = \'unknown\' )\n')
                        f.write(   '        do iz = 1, grd%nz\n')
                        f.write(   '          read(11,*) tmp, var%ni(isp,iz)\n')
                        f.write(   '        end do\n')
                        f.write(   '      close(11)\n')
                        f.write(   '    end if\n')

            f.write(   '\n')
            f.write(   '    var%ni_0 = var%ni\n')

            # boundary condition
            f.write(   '\n')
            f.write(   '    ! Lower boundary condition\n')
            f.write(   '    var%LowerBC = 0.0_dp\n')
            f.write(   '\n')
            for isp in range(len(bc_species)):
                if bc_Llabel[isp] != '4':
                    f.write(   "    isp = sp_index(spl, '"+bc_species[isp]+"')\n")
                    f.write(   '    var%LowerBC(isp,1) = '+str(bc_Llabel[isp])+'.0_dp\n')
                    f.write(   '    var%LowerBC(isp,2) = '+str(bc_Lval[isp])+'_dp\n')
                    f.write(   '\n')
            f.write(   '\n')
            f.write(   '    ! Upper boundary condition\n')
            f.write(   '    var%UpperBC = 0.0_dp\n')
            f.write(   '\n')
            for isp in range(len(bc_species)):
                if bc_Ulabel[isp] != '4':
                    f.write(   "    isp = sp_index(spl, '"+bc_species[isp]+"')\n")
                    if bc_Uval[isp] != 'jeans' and bc_Uval[isp] != 'Jeans':
                        f.write(   '    var%UpperBC(isp,1) = '+str(bc_Ulabel[isp])+'.0_dp\n')
                        f.write(   '    var%UpperBC(isp,2) = '+str(bc_Uval[isp])+'_dp\n')
                    if bc_Uval[isp] == 'jeans' or bc_Uval[isp] == 'Jeans':
                        f.write(   '    var%UpperBC(isp,1) = 10.0_dp ! Jeans escape\n')
                    f.write(   '\n')

            f.write(   '\n')
            f.write(   '  end subroutine v__in__exe\n')

            f.write(   '\n')
            f.write(   'end module v__in\n')

        print('Output f90 module "v__in.f90" in the directory '+Planet+'/'+dir0+' !')

        path = './CMakeLists.txt'
        with open(path, mode = 'w') as f:
            f.write('cmake_minimum_required(VERSION 3.19)\n')
            f.write('enable_language(Fortran) \n')
            f.write('project(PhotoChemistry)\n')
            f.write('set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)\n')
            f.write('add_executable(e__main\n')
            f.write('    e__main.f90\n')
            f.write('    v__tdec.f90\n')
            f.write('    c__prm.f90\n')
            f.write('    p__io.f90\n')
            f.write('    p__search.f90\n')
            f.write('    v__Venus.f90\n')
            f.write('    v__Earth.f90\n')
            f.write('    v__Mars.f90\n')
            f.write('    v__Jupiter.f90\n')
            f.write('    p__EUVAC.f90\n')
            f.write('    p__UV.f90\n')
            f.write('    p__eddy_diffusion.f90\n')
            f.write('    p__molecular_diffusion.f90\n')
            f.write('    p__photochem_opticaldepth.f90\n')
            f.write('    p__photochem_rate.f90\n')
            f.write('    p__photochem_transport.f90\n')
            f.write('    p__photochem_scheme.f90\n')
            f.write('    '+Planet+'/'+dir0+'/v__in.f90\n')
            f.write(')\n')
            f.write('\n')
            f.write('#ifort\n')
            f.write('if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)\n')
            #f.write('    set(CMAKE_Fortran_FLAGS         "-O3 -p")\n')
            #f.write('    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -warn all -standf95")\n')
            #f.write('    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")\n')
            f.write('endif()\n')
            f.write('#gfortran\n')
            f.write('if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)\n')
            #f.write('    set(CMAKE_Fortran_FLAGS         "-O3 -p -fbounds-check")\n')
            #f.write('    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -Wall")\n')
            #f.write('    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")\n')
            f.write('endif()\n')

        print('CMakeList updated!')
        
        path = './'+Planet+'/'+dir0+'/settings/plt_species.dat'
        with open(path, mode = 'w') as f:
            for isp in range(len(species)-1):
                f.write(species[isp]+'\n')
            f.write(species[len(species)-1])

        ###########################################################
        #
        #                Run PhotoChemical Model
        #
        ###########################################################
        if action == 'Run':
            path = './CMakeLists.txt'
            with open(path, mode = 'w') as f:
                f.write('cmake_minimum_required(VERSION 3.19)\n')
                f.write('enable_language(Fortran) \n')
                f.write('project(PhotoChemistry)\n')
                f.write('set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)\n')
                f.write('add_executable(e__main\n')
                f.write('    e__main.f90\n')
                f.write('    v__tdec.f90\n')
                f.write('    c__prm.f90\n')
                f.write('    p__io.f90\n')
                f.write('    p__search.f90\n')
                f.write('    v__Venus.f90\n')
                f.write('    v__Earth.f90\n')
                f.write('    v__Mars.f90\n')
                f.write('    v__Jupiter.f90\n')
                f.write('    p__EUVAC.f90\n')
                f.write('    p__UV.f90\n')
                f.write('    p__eddy_diffusion.f90\n')
                f.write('    p__molecular_diffusion.f90\n')
                f.write('    p__photochem_opticaldepth.f90\n')
                f.write('    p__photochem_rate.f90\n')
                f.write('    p__photochem_transport.f90\n')
                f.write('    p__photochem_scheme.f90\n')
                f.write('    '+Planet+'/'+dir0+'/v__in.f90\n')
                f.write(')\n')
                f.write('\n')
                f.write('#ifort\n')
                f.write('if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)\n')
                #f.write('    set(CMAKE_Fortran_FLAGS         "-O3 -p")\n')
                #f.write('    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -warn all -standf95")\n')
                #f.write('    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")\n')
                f.write('endif()\n')
                f.write('#gfortran\n')
                f.write('if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)\n')
                #f.write('    set(CMAKE_Fortran_FLAGS         "-O3 -p -fbounds-check")\n')
                #f.write('    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -Wall")\n')
                #f.write('    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")\n')
                f.write('endif()\n')

            print('Output f90 module "v__in.f90" in the directory '+Planet+'/'+dir0+' !')
            print('CMakeList updated!')

            # Run Model
            #os.system('cd build')
            #os.system('cmake ..')
            #os.system('cmake --build .')
            #os.system('cd ..')
            #os.system('./e__main')

            path = './e__main'
            if os.path.exists(path) == True:
                os.system('rm e__main')
            os.system('sh ./compile.sh')

            path = './'+Planet+'/'+dir0+'/settings/plt_species.dat'
            with open(path, mode = 'w') as f:
                for isp in range(len(species)-1):
                    f.write(species[isp]+'\n')
                f.write(species[len(species)-1])


# Make main Window ########################################################################################################
win = tk.Tk()
win.title("PROTEUS")
win.geometry("1280x800")
# This is the only one global variable in this python code
global reaction_win_label
reaction_win_label = 0

#checkbox(button) object needs a variable of BooleanVar
#this list stores a list of booleanVar(true or false) of each reaction
reaction_chk_bln = {} #koyama why is it dictionary?
for ich in range(len(reaction_rate_list)):
    reaction_chk_bln[ich] = tk.BooleanVar()
    reaction_chk_bln[ich].set(False)

#Set fix species window
#species is fixed or not
fix_species_bln = {}
input_species_char = {} #all species is input
for i in range(1000): #just take a large enough number
    fix_species_bln[i] = tk.BooleanVar()
    fix_species_bln[i].set(True)

search_species = {} #ex) ["CO2", "O"]
search_list = [1 for ich in range(len(reaction_rate_list))] #a list consisting of 1
#1 -> search
#0 -> not search

list_s = 0 #starting index number of reaction list
list_e = 0 #ending index number of reaction list

### callback functions ###
def callback_directory_window(iplnt,  list_s, list_e, reaction_chk_bln, fix_species_bln, input_species_char, search_list):
    def dummy():
        #Planet_list = [planet_name, starting index number of reaction list, ending index number of reaction list]
        Planet = Planet_list[iplnt][0] #planet_name
        list_s = Planet_list[iplnt][1] #starting index number of reaction list
        list_e = Planet_list[iplnt+1][1]-1 #ending index number of reaction list

        path = './'+Planet
        if os.path.exists(path) == False:
            os.makedirs(path)

        directory_window(iplnt, Planet, list_s, list_e, reaction_chk_bln, fix_species_bln, input_species_char, search_list)

    return dummy

#update reaction button (show all reaction button)
def callback_update_reaction_window(main_canvas, xbar, ybar, frame, upper_canvas, lower_canvas,
                                    iplnt, Planet, list_s, list_e, dir0, version,
                                    reaction_chk_bln, fix_species_bln, input_species_char,
                                    search_list):
    def dummy():
        main_canvas.destroy()
        xbar.destroy()
        ybar.destroy()
        frame.destroy()
        upper_canvas.destroy()
        lower_canvas.destroy()
        # return all the elements of search_list back to 1
        for ich in range(len(reaction_rate_list)):
            search_list[ich] = 1
        reaction_window(iplnt, Planet, list_s, list_e, dir0, version,
                        reaction_chk_bln, fix_species_bln, input_species_char,
                        search_list)
    return dummy

# search button
# this searching rule is not user-friendly koyama
def callback_search_reaction_list(main_canvas, xbar, ybar, frame, upper_canvas, lower_canvas,
                                  iplnt, Planet, list_s, list_e, dir0, version,
                                  reaction_chk_bln, fix_species_bln, input_species_char,
                                  search_species, search_text, search_list):
    def dummy():
        if search_text.get() != '':
            # AND search
            #seach_text = textbox that users put words in
            search_species = search_text.get().split('and') #it splits text by and
            if len(search_species) >= 2:
                for ich in range(len(reaction_rate_list)):
                    search_list[ich] = 1
                for isp in range(len(search_species)):
                    search_species[isp] = search_species[isp].lstrip()
                    search_species[isp] = search_species[isp].rstrip()
                    search_species[isp] = ' '+search_species[isp]+' '
                    for ich in range(len(reaction_rate_list)):
                        if ich >= list_s and ich <= list_e:
                            if search_species[isp] not in reaction_rate_list[ich]:
                                search_list[ich] = 0
            # OR search
            if len(search_species) == 1:
                search_species = search_text.get().split('or')
                for ich in range(len(reaction_rate_list)):
                    search_list[ich] = 0
                for isp in range(len(search_species)):
                    search_species[isp] = search_species[isp].lstrip()
                    search_species[isp] = search_species[isp].rstrip()
                    search_species[isp] = ' '+search_species[isp]+' '
                    for ich in range(len(reaction_rate_list)):
                        if ich >= list_s and ich <= list_e:
                            if search_species[isp] in reaction_rate_list[ich]:
                                search_list[ich] = 1

            main_canvas.destroy()
            xbar.destroy()
            ybar.destroy()
            frame.destroy()
            lower_canvas.destroy()
            upper_canvas.destroy()

            reaction_window(iplnt, Planet, list_s, list_e, dir0, version,
                            reaction_chk_bln, fix_species_bln, input_species_char,
                            search_list)
    return dummy

#reaction analysis
def callback_reaction_analysis(action, iplnt, reaction_chk_bln, fix_species_bln, dir0):
    def dummy():
        reaction_analysis(action, iplnt, reaction_chk_bln, fix_species_bln, dir0)
    return dummy

#Grid set button
def callback_grid_window(Planet, dir0):
    def dummy():
        dz_text = {}
        zrange_text = {}
        path = './'+Planet+'/'+dir0+'/settings/zgrid.dat'
        nz = 0
        if os.path.exists(path) == True:
            with open(path, mode = 'r') as f:
                lines = f.readlines()
            f.close()
            nz = len(lines)
        if nz == 0:
            nz = nz + 1
        Grid_win = tk.Toplevel()
        Grid_win.title("Set grid")
        Grid_win.geometry("600x800")
        grid_window(Planet, dir0, dz_text, zrange_text, nz, Grid_win)
    return dummy

def callback_add_grid_window(Planet, dir0, dz_text, zrange_text, nz, Grid_win, Grid_canvas):
    def dummy():
        path = './'+Planet+'/'+dir0+'/settings/zgrid.dat'
        with open(path, mode = 'w') as f:
            for iz in range(len(dz_text)):
                f.write(dz_text[iz].get()+', '+zrange_text[iz].get()+'\n')
        f.close()
        nz1 = nz + 1
        Grid_canvas.destroy()
        grid_window(Planet, dir0, dz_text, zrange_text, nz1, Grid_win)
    return dummy

def callback_remove_grid_window(Planet, dir0, dz_text, zrange_text, nz, Grid_win, Grid_canvas):
    def dummy():
        path = './'+Planet+'/'+dir0+'/settings/zgrid.dat'
        with open(path, mode = 'w') as f:
            for iz in range(len(dz_text)-1):
                f.write(dz_text[iz].get()+', '+zrange_text[iz].get()+'\n')
        f.close()
        nz1 = nz - 1
        Grid_canvas.destroy()
        grid_window(Planet, dir0, dz_text, zrange_text, nz1, Grid_win)
    return dummy

def callback_done_grid_window(Planet, dir0, dz_text, zrange_text, nz, Grid_win, Grid_canvas):
    def dummy():
        # read grid
        path = './'+Planet+'/'+dir0+'/settings/zgrid.dat'
        with open(path, mode = 'w') as f:
            for iz in range(nz):
                f.write(dz_text[iz].get()+', '+zrange_text[iz].get()+'\n')
        f.close()
        Grid_canvas.destroy()
        Grid_win.destroy()

    return dummy

def callback_input_window(Planet, dir0, fix_species_bln, input_species_char):
    def dummy():
        input_window(Planet, dir0, fix_species_bln, input_species_char)
    return dummy

def callback_enter_input_window(Planet, dir0, fix_species_bln, input_species_char,
                                input_win, input_list):
    def dummy():
        input_species = input_list.get().split(',')
        for i in range(len(input_species)):
            input_species[i] = input_species[i].lstrip().rstrip()
        input_win.destroy()
        input_fix_window(Planet, dir0, fix_species_bln, input_species_char,
                         search_list, input_species)
    return dummy

def callback_done_input_fix_window(Planet, dir0, fix_species_bln, input_species_char,
                                   dir_input_species, fname_input_species, input_species,
                                   input_fix_win, fix_canvas):
    def dummy():
        dir1 = dir_input_species.get().strip('\n')
        path = {}
        fname = {}
        for i in range(len(fname_input_species)):
            fname[i] = fname_input_species[i].get().strip('\n')
            path[i] = dir0+fname[i]
        print(fname)

        insp_all = ''
        for i in range(len(input_species)):
            if i < len(input_species)-1:
                insp_all += input_species[i]+', '
            if i == len(input_species)-1:
                insp_all += input_species[i]
        insp_all.rstrip()

        label = {}
        line = []
        path = './'+Planet+'/'+dir0+'/settings/input_path_list.dat'
        if os.path.exists(path) == True:
            with open(path, mode = 'r') as f:
                tmp = f.readline()
                tmp = f.readline()
                line = f.readlines()
            f.close()
        with open(path, mode = 'w') as f:
            f.write(insp_all+'\n')
            f.write(dir1+'\n')
            for j in range(len(input_species)):
                if fix_species_bln[j].get() == True:
                    bln = '1'
                if fix_species_bln[j].get() == False:
                    bln = '0'
                f.write(input_species[j]+':'+fname[j]+';'+bln+'\n')
            for i in range(len(line)):
                sp = re.findall('(.*):',line[i])
                label[i]=0
                for j in range(len(input_species)):
                    if input_species[j] == sp[0].lstrip().rstrip():
                        label[i]=1
            for i in range(len(line)):
                if label[i] == 0:
                    f.write(line[i])
        f.close()


        fix_canvas.destroy()
        input_fix_win.destroy()
    return dummy

def callback_help_window():
    def dummy():
        help_window()
    return dummy

def callback_exit_window(win):
    def dummy():
        win.destroy()
    return dummy

def callback_create_dir_window(iplnt, Planet, list_s, list_e,
                               reaction_chk_bln, fix_species_bln, input_species_char, search_list,
                               dir_win, dir_create):
    def dummy():

        dir0 = dir_create.get().strip('\n').strip('/')
        dir_win.destroy()

        if dir0 == '':
            err_win = tk.Toplevel() #display on the main window
            err_win.title("error")
            err_win.geometry("300x100")
            txt='Please enter the name of directory again!'
            txt_label = tk.Label(err_win, anchor="w", justify="left", font=('', 15), text=txt)
            txt_label.place(x=10, y=10)

        if dir0 != '':

            path = './'+Planet+'/dir_list.dat'
            if os.path.exists(path) == True:
                with open(path, mode = 'r') as f:
                    dir_list = f.readlines()
                f.close()
                with open(path, mode = 'w') as f:
                    for i in range(len(dir_list)):
                        f.write(dir_list[i])
                    f.write('\n')
                    f.write(dir0)
                f.close()
            if os.path.exists(path) == False:
                with open(path, mode = 'w') as f:
                    f.write(dir0)
                f.close()

            path = './'+Planet
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/input/PLJ_list'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/input/density'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/input/Temperature'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/settings'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/info'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/output'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/output/density'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/output/density/3Drot'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/output/density/2Dstable'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/output/density/num'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/output/density/vmr'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/output/flux'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir0+'/output/rate'
            if os.path.exists(path) == False:
                os.makedirs(path)

            reaction_window(iplnt, Planet, list_s, list_e, dir0, version, reaction_chk_bln, fix_species_bln, input_species_char, search_list)

    return dummy

def callback_select_dir_window(iplnt, Planet, list_s, list_e, dir0, version, reaction_chk_bln, fix_species_bln, input_species_char, search_list, dir_win):
    def dummy():

        dir1 = dir0.strip('\n')

        path = './'+Planet+'/'+dir1
        if os.path.exists(path) ==True:
            dir_win.destroy()

            path = './'+Planet+'/'+dir1+'/input/PLJ_list'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/input/density'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/input/Temperature'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/settings'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/info'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/output'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/output/density'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/output/density/3Drot'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/output/density/3Drot/all'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/output/density/2Dstable'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/output/density/num'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/output/density/vmr'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/output/flux'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir1+'/output/rate'
            if os.path.exists(path) == False:
                os.makedirs(path)

            reaction_window(iplnt, Planet, list_s, list_e, dir1, version, reaction_chk_bln, fix_species_bln, input_species_char, search_list)

        path = './'+Planet+'/'+dir1
        if os.path.exists(path) ==False:
            dir_win.destroy()
            err_win = tk.Toplevel() #display on the main window
            err_win.title("error")
            err_win.geometry("300x100")
            txt='Directory "'+dir1+'" does not exist!\nPlease select again!'
            txt_label = tk.Label(err_win, anchor="w", justify="left", font=('', 15), text=txt)
            txt_label.place(x=10, y=10)

            path1 = './'+Planet+'/dir_list.dat'
            if os.path.exists(path1) == True:
                with open(path1, mode = 'r') as f:
                    dir_list = f.readlines()
                f.close()
                with open(path1, mode = 'w') as f:
                    for i in range(len(dir_list)):
                        dir_list[i] = dir_list[i].strip('\n')
                        if dir_list[i] != dir1:
                            f.write(dir_list[i])
                            if i < len(dir_list)-2:
                                f.write('\n')
                f.close()
            directory_window(iplnt, Planet, list_s, list_e, reaction_chk_bln, fix_species_bln, input_species_char, search_list)

    return dummy

def callback_save_dir_window(iplnt, Planet, list_s, list_e, dir0, version,
                             reaction_chk_bln, fix_species_bln, input_species_char, search_list,
                             dir_win, dir_save):
    def dummy():

        dir1 = dir0.strip('\n')
        dir2 = dir_save.get().strip('\n')

        path = './'+Planet+'/'+dir1
        if os.path.exists(path) ==True:
            dir_win.destroy()
            os.system('cp -R ./'+Planet+'/'+dir1+' ./'+Planet+'/'+dir2)

            path = './'+Planet+'/dir_list.dat'
            if os.path.exists(path) == True:
                with open(path, mode = 'r') as f:
                    dir_list = f.readlines()
                f.close()
                with open(path, mode = 'w') as f:
                    for i in range(len(dir_list)):
                        f.write(dir_list[i])
                    f.write('\n')
                    f.write(dir2)
                f.close()
            if os.path.exists(path) == False:
                with open(path, mode = 'w') as f:
                    f.write(dir2)
                f.close()

            path = './'+Planet+'/'+dir2+'/input/PLJ_list'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir2+'/input/density'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir2+'/input/Temperature'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir2+'/settings'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir2+'/info'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir2+'/output'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir2+'/output/density'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir2+'/output/density/3Drot'
            if os.path.exists(path) == False:
                os.makedirs(path)
            
            path = './'+Planet+'/'+dir2+'/output/density/2Dstable'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir2+'/output/density/num'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir2+'/output/density/vmr'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir2+'/output/flux'
            if os.path.exists(path) == False:
                os.makedirs(path)

            path = './'+Planet+'/'+dir2+'/output/rate'
            if os.path.exists(path) == False:
                os.makedirs(path)

            # rename directory path of fixed species
            path = './'+Planet+'/'+dir2+'/settings/input_path_list.dat'
            if os.path.exists(path) == True:
                with open(path, mode = 'r') as f:
                    lines = f.readlines()
                f.close()
                with open(path, mode = 'w') as f:
                    f.write(lines[0])
                    f.write(dir2+'/input/density\n')
                    for i in range(len(lines)):
                        if i >= 2:
                            f.write(lines[i])
                f.close()

            reaction_window(iplnt, Planet, list_s, list_e, dir2, version, reaction_chk_bln, fix_species_bln, input_species_char, search_list)

        path = './'+Planet+'/'+dir1
        if os.path.exists(path) ==False:
            dir_win.destroy()
            err_win = tk.Toplevel() #display on the main window
            err_win.title("error")
            err_win.geometry("300x100")
            txt='Directory "'+dir1+'" does not exist!\nPlease try again!'
            txt_label = tk.Label(err_win, anchor="w", justify="left", font=('', 15), text=txt)
            txt_label.place(x=10, y=10)

            path1 = './'+Planet+'/dir_list.dat'
            if os.path.exists(path1) == True:
                with open(path1, mode = 'r') as f:
                    dir_list = f.readlines()
                f.close()
                with open(path1, mode = 'w') as f:
                    for i in range(len(dir_list)):
                        dir_list[i] = dir_list[i].strip('\n')
                        if dir_list[i] != dir1:
                            f.write(dir_list[i])
                            if i < len(dir_list)-2:
                                f.write('\n')
                f.close()
            directory_window(iplnt, Planet, list_s, list_e, reaction_chk_bln, fix_species_bln, input_species_char, search_list)

    return dummy

def callback_boundary_condition_input_window(Planet, dir0):
    def dummy():
        boundary_condition_input_window(Planet, dir0)
    return dummy

def callback_enter_boundary_condition_input_window(Planet, dir0, bc_input_win, bc_list):
    def dummy():
        bc_species = bc_list.get().split(',')
        for i in range(len(bc_species)):
            bc_species[i] = bc_species[i].lstrip().rstrip()
        bc_input_win.destroy()
        boundary_condition_set_window(Planet, dir0, bc_species)
    return dummy

def callback_done_boundary_condition_set_window(Planet, dir0, bc_species, bcvl, bcvu, bcset_win, 
                                                lbc_val1, lbc_val2, lbc_val3,
                                                ubc_val1, ubc_val2, ubc_val3):
    def dummy():

        path = './'+Planet+'/'+dir0+'/settings/boundary_condition.dat'
        if os.path.exists(path) == True:
            with open(path, mode = 'r') as f:
                tmp = f.readline()
                line = f.readlines()
            f.close()
        with open(path, mode = 'w') as f:
            for i in range(len(bc_species)-1):
                f.write(bc_species[i]+',')
            f.write(bc_species[len(bc_species)-1])
            f.write('\n')
            for i in range(len(bc_species)):
                nl = bcvl[i].get()//10000
                nu = bcvu[i].get()//10000
                f.write(bc_species[i]+' :')
                if nl == 1:
                    f.write(' L '+str(nl)+' lval '+lbc_val1[i].get())
                if nl == 2:
                    f.write(' L '+str(nl)+' lval '+lbc_val2[i].get())
                if nl == 3:
                    f.write(' L '+str(nl)+' lval '+lbc_val3[i].get())
                if nl == 4:
                    f.write(' L '+str(nl)+' lval '+'0')
                if nu == 1:
                    f.write(' U '+str(nu)+' uval '+ubc_val1[i].get())
                if nu == 2:
                    f.write(' U '+str(nu)+' uval '+ubc_val2[i].get())
                if nu == 3:
                    f.write(' U '+str(nu)+' uval '+ubc_val3[i].get())
                if nu == 4:
                    f.write(' U '+str(nu)+' uval '+'0')
                f.write('\n')

        bcset_win.destroy()

    return dummy

def callback_calculation_set_window(Planet, dir0):
    def dummy():
        calculation_set_window(Planet, dir0)
    return dummy

def callback_done_calculation_set_window(Planet, dir0, 
                                         rbval, 
                                         LAT, f107, ls, 
                                         stepmax, 
                                         sumdt, dtmax, 
                                         sumdtrot, dtmaxrot, 
                                         calcset_win):
    def dummy():

        path = './'+Planet+'/'+dir0+'/settings/calculation_setting.dat'
        if os.path.exists(path) == True:
            with open(path, mode = 'r') as f:
                line = f.readlines()
            f.close()
        with open(path, mode = 'w') as f:
            f.write('Dimension : '+ str(rbval[1].get())+'\n')
            f.write('LAT : '+ str(LAT.get())+'\n')
            f.write('F10.7 : '+ str(f107.get())+'\n')
            f.write('Ls : '+ str(ls.get())+'\n')
            f.write('stepmax : '+ str(stepmax.get())+'\n')
            f.write('sumdtst : '+ str(sumdt.get())+'\n')
            f.write('dtmaxst : '+ str(dtmax.get())+'\n')
            f.write('sumdtrot : '+ str(sumdtrot.get())+'\n')
            f.write('dtmaxrot : '+ str(dtmaxrot.get())+'\n')
            f.write('scheme : '+ str(rbval[2].get())+'\n')
            f.write('inversion : '+ str(rbval[3].get())+'\n')

        calcset_win.destroy()
    return dummy

def callback_plot_window(Planet, dir0):
    def dummy():
        plot_window(Planet, dir0)
    return dummy


def callback_plot(Planet, dir0, species, action, rbvar, adv, fs, 
                  sp_chk_bln, sp1, sp2, 
                  xr, yr, LATin2, LATin3, LTin3, 
                  savefname, save_bln, ftrans_bln):

    def dummy():
        path = './'+Planet+'/'+dir0+'/settings/plt_range.dat'
        with open(path, mode = 'w') as f:
            f.write('xrm3:'+xr[0][0].get()+':'+xr[0][1].get()+'\n')
            f.write('xrcm3:'+xr[1][0].get()+':'+xr[1][1].get()+'\n')
            f.write('xrvmr:'+xr[2][0].get()+':'+xr[2][1].get()+'\n')
            f.write('yr:'+yr[0].get()+':'+yr[1].get())

        path = './'+Planet+'/'+dir0+'/settings/plt_fontsize.dat'
        with open(path, mode = 'w') as f:
            f.write('xlabel:'+fs[0].get()+'\n')
            f.write('ylabel:'+fs[1].get()+'\n')
            f.write('tick:'+fs[2].get()+'\n')
            f.write('legend:'+fs[3].get()+'\n')
            f.write('linewidth:'+fs[4].get())

        if rbvar.get() == 1:
            action1 = action
            LAT = 0
            LT = 0

        if rbvar.get() == 2:
            action1 = '2D Lat '+action
            LAT = float(LATin2.get().lstrip().rstrip())
            LT = 0

        if rbvar.get() == 3:
            action1 = '3D Rot '+action
            LAT = float(LATin3.get().lstrip().rstrip())
            LT = float(LTin3.get().lstrip().rstrip())

        savelabel = 0
        if save_bln.get() == True:
            savelabel = 1
            if ftrans_bln.get() == True:
                savelabel = 2
            path = './'+Planet+'/'+dir0+'/figure'
            if os.path.exists(path) == False:
                os.makedirs(path)

        plot(Planet, dir0, species, action1, adv, fs, 
             sp_chk_bln, sp1, sp2, 
             xr, yr, LAT, LT, 
             savelabel, savefname)

    return dummy

def callback_plot_order(Planet, dir0):
    def dummy():
        plot_order_window(Planet, dir0)
    return dummy

def callback_plot_order_done(Planet, dir0, text, plt_win):
    def dummy():
        path = './'+Planet+'/'+dir0+'/settings/plt_species_order.dat'
        with open(path, mode = 'w') as f:
            f.write(text.get(1.0, tk.END))
        plt_win.destroy()
    return dummy

### Window definitions ###

# Directory window
def directory_window(iplnt, Planet, list_s, list_e, reaction_chk_bln, fix_species_bln, input_species_char, search_list):

    dir_win = tk.Toplevel() #display on the main window
    dir_win.title("Select Project")
    dir_win.geometry("720x600")

    dir_list = []
    path = './'+Planet+'/dir_list.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            dir_list = f.readlines()
        f.close()

    dir_canvas = tk.Canvas(dir_win, width=720,height=100, highlightthickness=0)

    ybar = tk.Scrollbar(dir_win, orient=tk.VERTICAL) #scroll bar
    ybar.pack(side=tk.RIGHT, fill=tk.Y)
    ybar.config(command=dir_canvas.yview)

    #frame on the main canvas
    dir_frame = tk.Frame(dir_canvas, width=720, height=500 + len(dir_list) * 30)
    dir_canvas.create_window((0,0), window=dir_frame, anchor=tk.NW, width=dir_canvas.cget('width')) #place frame on the canvas

    dir_canvas.config(yscrollcommand=ybar.set)
    dir_canvas.config(scrollregion=(0,0,720,300 + len(dir_list) * 30))
    dir_canvas.pack(anchor=tk.NW, expand=1, fill=tk.BOTH)
    dir_canvas.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    dir_frame.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    txt  = Planet
    txt_label = tk.Label(dir_frame, anchor="w", justify="left", font=('', 20), text=txt)
    txt_label.place(x=10, y=10)
    txt_label.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    txt  = '--------------------------\n'
    txt += 'Select or Create Project'
    txt_label = tk.Label(dir_frame, anchor="w", justify="left", font=('', 15), text=txt)
    txt_label.place(x=10, y=40)
    txt_label.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    txt  = '    # If you want to create new directory, please enter the name of new project (derectory) and press "Create New".\n'
    txt += '    # If you want to use a directory that already exists, please pless "select" to select from the project (directory) list below. \n'
    txt += '    # If you want to load a directory and save as new project (directory), \n'
    txt += '        please enter the name of new project (directory) name and press "Save as".\n'
    txt += '    # Location of the project (directory) is: ./'+Planet+'/"project (directory) name"\n'
    txt_label = tk.Label(dir_frame, anchor="w", justify="left", font=('', 12), text=txt)
    txt_label.place(x=10, y=90)
    txt_label.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    yline = 180
    dir_create = tk.Entry(dir_frame, width=20)
    dir_create.place(x=90, y=yline)
    dir_create.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    create_btn = tk.Button(dir_frame, font=('', 12), text=u'Create New')
    create_btn.place(x=10, y=yline+4)
    create_btn.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    create_btn["command"] = callback_create_dir_window(iplnt, Planet, list_s, list_e,
                                                       reaction_chk_bln, fix_species_bln, input_species_char, search_list,
                                                       dir_win, dir_create)

    txt  = 'Project (directory) list:'
    txt_label = tk.Label(dir_frame, anchor="w", justify="left", font=('', 12), text=txt)
    txt_label.place(x=10, y=yline+40)
    txt_label.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
 
    for i in range(len(dir_list)):
        dir0 = dir_list[i]

        if dir0.strip('\n') != '':
            select_btn = tk.Button(dir_frame, font=('', 12), text=u'Select')
            select_btn.place(x=10, y=yline+70+(i*30))
            select_btn.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
            select_btn["command"] = callback_select_dir_window(iplnt, Planet, list_s, list_e, dir0, version,
                                                                reaction_chk_bln, fix_species_bln, input_species_char, search_list,
                                                                dir_win)
            dir_label = tk.Label(dir_frame, font=('', 12), text=dir0)
            dir_label.place(x=70, y=yline+70+(i*30))
            dir_label.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

            dir_save = tk.Entry(dir_frame, width=20)
            dir_save.place(x=360, y=yline+66+(i*30))
            dir_save.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
            save_btn = tk.Button(dir_frame, font=('', 12), text=u'Save as')
            save_btn.place(x=300, y=yline+70+(i*30))
            save_btn.bind("<MouseWheel>", lambda e:dir_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
            save_btn["command"] = callback_save_dir_window(iplnt, Planet, list_s, list_e, dir0, version,
                                                            reaction_chk_bln, fix_species_bln, input_species_char, search_list,
                                                            dir_win, dir_save)


# Help window # Manual
def help_window():
    help_win = tk.Toplevel() #display on the main window
    help_win.title("Help information")
    help_win.geometry("1280x800")

    text = tk.Text(help_win, font=("",15), height=200, width=190, highlightthickness=0)
    text.pack()

    text.insert(tk.END, '############################################\n')
    text.insert(tk.END, ' Help information\n')
    text.insert(tk.END, '############################################\n')
    text.insert(tk.END, '\n')
    text.insert(tk.END, ' # Search window help\n')
    text.insert(tk.END, '   * AND search\n')
    text.insert(tk.END, '   -> species should be separated by \' and \'\n')
    text.insert(tk.END, '   (ex) \'CO2 and O\'      ,      \'CO and CO+\'\n')
    text.insert(tk.END, '\n')
    text.insert(tk.END, '   * OR search\n')
    text.insert(tk.END, '   -> species should be separated by \' or \'\n')
    text.insert(tk.END, '   (ex) \'CO2 or O\'      ,      \'CO or CO+\'\n')

# Grid window
def grid_window(Planet, dir0, dz_text, zrange_text, nz, Grid_win):

    set_dz = {}
    set_zrange = {}
    dz = {}
    zrange = {}

    Grid_canvas = tk.Canvas(Grid_win, width=600,height=800)
    Grid_canvas.create_rectangle(0, 0, 600, 800)
    Grid_canvas.place(x=0, y=0)

    path = './'+Planet+'/'+dir0+'/settings/zgrid.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            lines = f.readlines()
            for iz in range(len(lines)):
                lchar = lines[iz]
                rchar = lines[iz]
                lchar = lchar.split(',')
                dz[iz] = lchar[0]
                rchar = rchar.split(',')
                rchar = rchar[1].split('\n')
                zrange[iz] = rchar[0].lstrip()
        f.close()

    for iz in range(nz):

        set_dz[iz] = tk.Label(Grid_canvas, anchor="w", text = u"dz [km]:")
        set_dz[iz].pack()
        set_dz[iz].place(x = 10, y = 83 + iz * 30 )
        dz_text[iz] = tk.Entry(Grid_win, width=5)
        if iz < len(dz):
            dz_text[iz].insert(tk.END, dz[iz])
        dz_text[iz].pack()
        dz_text[iz].place(x=70,y = 80 + iz * 30 )

        set_zrange[iz] = tk.Label(Grid_canvas, anchor="w" ,text = u"z range [km]:")
        set_zrange[iz].pack()
        set_zrange[iz].place(x = 150, y = 83 + iz * 30 )
        zrange_text[iz] = tk.Entry(Grid_win, width=20)
        if iz < len(dz):
            zrange_text[iz].insert(tk.END, zrange[iz])
        zrange_text[iz].pack()
        zrange_text[iz].place(x=250,y = 80 + iz * 30 )

    Add_Grid_btn = tk.Button(Grid_canvas, text=u'Add')
    Add_Grid_btn["command"] = callback_add_grid_window(Planet, dir0, dz_text, zrange_text, nz, Grid_win, Grid_canvas)
    Add_Grid_btn.pack()
    Add_Grid_btn.place(x=450, y = 83 + nz * 30 - 30)

    Remv_Grid_btn = tk.Button(Grid_canvas, text=u'Remove')
    Remv_Grid_btn["command"] = callback_remove_grid_window(Planet, dir0, dz_text, zrange_text, nz, Grid_win, Grid_canvas)
    Remv_Grid_btn.pack()
    Remv_Grid_btn.place(x=490, y = 83 + nz * 30 - 30)

    Done_btn = tk.Button(Grid_canvas, text=u'Done')
    Done_btn["command"] = callback_done_grid_window(Planet, dir0, dz_text, zrange_text, nz, Grid_win, Grid_canvas)
    Done_btn.pack()
    Done_btn.place(x=550, y = 83 + nz * 30 - 30)

# Set input species window
def input_window(Planet, dir0, fix_species_bln, input_species_char):
    input_win = tk.Toplevel()
    input_win.title("Set input species")
    input_win.geometry("400x200")
    char0 = tk.Label(input_win, font=('',20), anchor="w", text = u"Set input species")
    char0.pack()
    char0.place(x=10, y=10)
    char1 = tk.Label(input_win, anchor="w", justify="left", text = u"Please enter input species below.\nEach species must be devided by ','.")
    char1.pack()
    char1.place(x=10, y=40)

    input_list = tk.Entry(input_win, width=30)
    input_list.place(x=10, y=100)

    # read path of fixed species
    path = './'+Planet+'/'+dir0+'/settings/input_path_list.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            line = f.readline()
            line = line.strip('\n')
            input_list.insert(tk.END, line)
        f.close()
    if os.path.exists(path) == False:
        with open(path, mode = 'w') as f:
            f.write('')
        f.close()

    enter_btn = tk.Button(input_win, font=('', 12), text=u'enter')
    enter_btn.place(x=300, y=105)
    enter_btn["command"] = callback_enter_input_window(Planet, dir0,fix_species_bln, input_species_char,
                                                       input_win, input_list)

# Set input path and Select fix window
def input_fix_window(Planet, dir0, fix_species_bln, input_species_char,
                     search_list, input_species):

    yline = 250

    input_fix_win = tk.Toplevel()
    input_fix_win.title("Set input species")
    input_fix_win.geometry("600x800")
    fix_canvas = tk.Canvas(input_fix_win, width=600,height=100, highlightthickness=0)

    ybar = tk.Scrollbar(input_fix_win, orient=tk.VERTICAL) #scroll bar
    ybar.pack(side=tk.RIGHT, fill=tk.Y)
    ybar.config(command=fix_canvas.yview)

    fix_canvas.config(yscrollcommand=ybar.set)
    fix_canvas.config(scrollregion=(0,0,600,yline + 120 + len(input_species) * 60))
    fix_canvas.pack(anchor=tk.NW, expand=1, fill=tk.BOTH)

    #frame on the main canvas
    fix_frame = tk.Frame(fix_canvas, width=600, height=yline + 120 + len(input_species) * 60)
    fix_canvas.create_window((0,0), window=fix_frame, anchor=tk.NW, width=fix_canvas.cget('width')) #place frame on the canvas
    fix_canvas.bind("<MouseWheel>", lambda e:fix_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    fix_frame.bind("<MouseWheel>", lambda e:fix_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    char0 = tk.Label(fix_frame, anchor="w", font=('', 20), text = u"Set path of input density files")
    char0.place(x = 10, y = 10 )
    char0.bind("<MouseWheel>", lambda e:fix_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char1 = tk.Label(fix_frame, anchor="w", justify='left', text = u"Please enter the directory of input density files below.\nDirectory is loacted at ./"+Planet+"/'directory name'")
    char1.place(x = 10, y = 40 )
    char1.bind("<MouseWheel>", lambda e:fix_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    dir_input_species = tk.Entry(fix_frame, width=30)
    dir_input_species.place(x=10, y=100)
    dir_input_species.bind("<MouseWheel>", lambda e:fix_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    # display
    char2 = tk.Label(fix_frame, anchor="w", text = u"Please check below if you need to fix the density:") #desity -> density koyama
    char2.place(x = 10, y = yline-80 )
    char2.bind("<MouseWheel>", lambda e:fix_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    # make fixed species check button
    chk_fix = {}
    fname_input_species = {}
    path_read = {}
    # read path of fixed species
    path = './'+Planet+'/'+dir0+'/settings/input_path_list.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            line = f.readline()
            dir1 = f.readline()
            dir_input_species.insert(tk.INSERT, dir1.strip('\n'))
            path_read = f.readlines()
            for isp in range(len(path_read)):
                path_read[isp] = path_read[isp].strip('\n')
        f.close()

        if dir1 == '':
            dir1 = './'+dir0+'/input/density'
            dir_input_species.insert(tk.INSERT, dir1)

    if os.path.exists(path) == False:
        dir1 = './'+dir0+'/input/density'
        dir_input_species.insert(tk.INSERT, dir1)

    for isp in range(len(input_species)):
        input_species_unicode = reaction_unicode(input_species[isp])
        chk_fix[isp] = tk.Checkbutton(fix_frame, width = 20, anchor="w", variable = fix_species_bln[isp], text = 'Fix '+input_species_unicode+' density')
        chk_fix[isp].place(x = 10, y = yline + isp * 60 )
        chk_fix[isp].bind("<MouseWheel>", lambda e:fix_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
        path = tk.Label(fix_frame, anchor="w", text = u"path:")
        path.place(x = 10, y = yline + 30 + isp * 60 )
        path.bind("<MouseWheel>", lambda e:fix_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        fname1 = ''
        fixlabel = '0'
        for j in range(len(path_read)):
            sp = re.findall('(.*):', path_read[j])
            fname0 = re.findall(':(.*)', path_read[j])
            if ';' in fname0[0]:
                fixlabel = re.findall(';(.*)', fname0[0])
                fixlabel = fixlabel[0].lstrip().rstrip()
                fname0 = re.findall('(.*);', fname0[0])
            if sp[0] == input_species[isp]:
                fname1 = fname0[0]
                fix_species_bln[isp].set(False)
                if fixlabel == '1':
                    fix_species_bln[isp].set(True)
                if fixlabel == '0':
                    fix_species_bln[isp].set(False)

        if fname1 == '':
            fname1 = input_species[isp]+'.dat'

        fname_input_species[isp] = tk.Entry(fix_frame, width=40)
        fname_input_species[isp].insert(tk.END, fname1)
        fname_input_species[isp].place(x=60,y = yline + 30 + isp * 60 )
        fname_input_species[isp].bind("<MouseWheel>", lambda e:fix_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    Done_btn = tk.Button(input_fix_win, font=('',15), text=u'Done')
    Done_btn["command"] = callback_done_input_fix_window(Planet, dir0, fix_species_bln, input_species_char,
                                                         dir_input_species, fname_input_species, input_species,
                                                         input_fix_win, fix_canvas)
    Done_btn.place(x=480, y = 10)

    # set all the check buttons in displayed window as True
    def all_Select_click():
        for isp in range(len(input_species)):
            fix_species_bln[isp].set(True)
    all_Select_Button = tk.Button(fix_frame, text='All Fix',command = all_Select_click)
    all_Select_Button.place(x=10, y=210)

    # set all the check buttons in displayed window as False
    def all_Clear_click():
        for isp in range(len(input_species)):
            fix_species_bln[isp].set(False)
    all_Clear_Button = tk.Button(fix_frame, text='All Variable',command = all_Clear_click)
    all_Clear_Button.place(x=60, y=210)

# Set boundary condition window
def boundary_condition_input_window(Planet, dir0):
    bc_input_win = tk.Toplevel()
    bc_input_win.title("Entry species for setting boundary condition")
    bc_input_win.geometry("400x200")
    char0 = tk.Label(bc_input_win, font=('',20), anchor="w", text = u"Entry species for setting boundary condition")
    char0.place(x=10, y=10)
    char1 = tk.Label(bc_input_win, anchor="w", justify="left", text = u"Please enter species below.\nEach species must be devided by ','.")
    char1.place(x=10, y=40)

    bc_list = tk.Entry(bc_input_win, width=30)
    bc_list.place(x=10, y=100)

    # read boundary conditions
    path = './'+Planet+'/'+dir0+'/settings/boundary_condition.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            line = f.readline()
            line = line.strip('\n')
            bc_list.insert(tk.END, line)
        f.close()
    if os.path.exists(path) == False:
        with open(path, mode = 'w') as f:
            f.write('')
        f.close()

    enter_btn = tk.Button(bc_input_win, font=('', 12), text=u'enter')
    enter_btn.place(x=300, y=105)
    enter_btn["command"] = callback_enter_boundary_condition_input_window(Planet, dir0, bc_input_win, bc_list)

# Set input path and Select fix window
def boundary_condition_set_window(Planet, dir0, bc_species):

    bcset_win = tk.Toplevel()
    bcset_win.title("Set boundary condition")
    bcset_win.geometry("720x800")

    # make type of boundary condition check button
    bc_read = {}
    bc_species_read = {}
    Llabel  = {}
    Lval    = {}
    Ulabel  = {}
    Uval    = {}
    for isp in range(len(bc_species)):
        Llabel[isp] = '0'
        Ulabel[isp] = '0'
        Lval[isp] = '0'
        Uval[isp] = '0'
    # read path of fixed species
    path = './'+Planet+'/'+dir0+'/settings/boundary_condition.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            line = f.readline()
            bc_read = f.readlines()
            for isp in range(len(bc_read)):
                bc_read[isp] = bc_read[isp].strip('\n')
                bcchar = re.findall(':(.*)',bc_read[isp])
                bc_species_read[isp]   = re.findall('(.*):',bc_read[isp])[0].lstrip().rstrip()
                Llabel[isp] = re.findall('L(.*)lval',bcchar[0])[0].lstrip().rstrip()
                Lval[isp]   = re.findall('lval(.*)U',bcchar[0])[0].lstrip().rstrip()
                Ulabel[isp] = re.findall('U(.*)uval',bcchar[0])[0].lstrip().rstrip()
                Uval[isp]   = re.findall('uval(.*)',bcchar[0])[0].lstrip().rstrip()
        f.close()
        

    bc_canvas = tk.Canvas(bcset_win, width=720,height=800,highlightthickness=0)

    ybar = tk.Scrollbar(bcset_win, orient=tk.VERTICAL) #scroll bar
    ybar.pack(side=tk.RIGHT, fill=tk.Y)
    ybar.config(command=bc_canvas.yview)

    #frame on the main canvas
    bc_frame = tk.Frame(bc_canvas, width=720, height=300 + len(bc_species) * 150)
    bc_canvas.create_window((0,0), window=bc_frame, anchor=tk.NW, width=bc_canvas.cget('width')) #place frame on the canvas

    text = tk.Text(bc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.pack()
    text.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '====================================\n')
    text.insert(tk.END, '  Boundary Condition Setting\n')
    text.insert(tk.END, '====================================\n')
    text.insert(tk.END, ' Please Check and insert values to set boundary condition.\n')
    text.insert(tk.END, '  # [n]: Boundary condition of density [/m\u00B3] \n')
    text.insert(tk.END, '  # [f]: Boundary condition of upward flux [/m\u00B2/s] \n')
    text.insert(tk.END, '  # [v]: Boundary condition of upward velocity [m/s] \n')

    bc_canvas.config(yscrollcommand=ybar.set)
    bc_canvas.config(scrollregion=(0,0,720,300 + len(bc_species) * 150))
    bc_canvas.pack(anchor=tk.NW, expand=1, fill=tk.BOTH)
    bc_canvas.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    bc_frame.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    ys = 200
    lower = tk.Label(bc_frame, anchor="w", justify='left', text = u"Lower Boundary")
    lower.place(x = 100, y = ys-30 )
    lower.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    upper = tk.Label(bc_frame, anchor="w", justify='left', text = u"Upper Boundary")
    upper.place(x = 350, y = ys-30 )
    upper.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    bcvl = {}
    bcvu = {}
    lbc_val1 = {}
    lbc_val2 = {}
    lbc_val3 = {}
    ubc_val1 = {}
    ubc_val2 = {}
    ubc_val3 = {}
    for i in range(len(bc_species)):
        jsp = -1
        for j in range(len(bc_species_read)):
            if bc_species[i] == bc_species_read[j]:
                jsp = j
        bc_species_unicode = reaction_unicode(bc_species[i])
        bcsp = tk.Label(bc_frame, anchor="w", justify='left', text = bc_species_unicode)
        bcsp.place(x = 40, y = ys + i*150 )
        bcsp.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        bcvl[i] = tk.IntVar()
        if jsp >= 0:
            bcvl[i].set(int(Llabel[jsp])*10000+i)
        if jsp < 0:
            bcvl[i].set(40000+i)

        rb1l = tk.Radiobutton(bc_frame,variable=bcvl[i],value=10000+i,text=u"[n]")
        rb1l.place(x=100, y = ys + i*150)
        rb1l.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
        rb2l = tk.Radiobutton(bc_frame,variable=bcvl[i],value=20000+i,text=u"[f]")
        rb2l.place(x=100, y = ys + i*150+30)
        rb2l.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
        rb3l = tk.Radiobutton(bc_frame,variable=bcvl[i],value=30000+i,text=u"[v]")
        rb3l.place(x=100, y = ys + i*150+60)
        rb3l.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
        rb4l = tk.Radiobutton(bc_frame,variable=bcvl[i],value=40000+i,text=u"non")
        rb4l.place(x=100, y = ys + i*150+90)
        rb4l.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        chl1 = tk.Label(bc_frame,text=u"[/m\u00B3]")
        chl1.place(x=280, y = ys + i*150)
        chl1.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
        chl2 = tk.Label(bc_frame,text=u"[/m\u00B2/s]")
        chl2.place(x=280, y = ys + i*150+30)
        chl2.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
        chl3 = tk.Label(bc_frame,text=u"[/m/s]")
        chl3.place(x=280, y = ys + i*150+60)
        chl3.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        lbc_val1[i] = tk.Entry(bc_frame, width=13)
        if jsp >= 0:
            if Llabel[jsp] == '1':
                lbc_val1[i].insert(tk.END, Lval[jsp])
        lbc_val1[i].place(x=145,y = ys + i*150 )
        lbc_val1[i].bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        lbc_val2[i] = tk.Entry(bc_frame, width=13)
        if jsp >= 0:
            if Llabel[jsp] == '2':
                lbc_val2[i].insert(tk.END, Lval[jsp])
        lbc_val2[i].place(x=145,y = ys + i*150+30 )
        lbc_val2[i].bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        lbc_val3[i] = tk.Entry(bc_frame, width=13)
        if jsp >= 0:
            if Llabel[jsp] == '3':
                lbc_val3[i].insert(tk.END, Lval[jsp])
        lbc_val3[i].place(x=145,y = ys + i*150+60 )
        lbc_val3[i].bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        
        bcvu[i] = tk.IntVar()
        if jsp >= 0:
            bcvu[i].set(int(Ulabel[jsp])*10000+i)
        if jsp < 0:
            bcvu[i].set(40000+i)

        rb1u = tk.Radiobutton(bc_frame,variable=bcvu[i],value=10000+i,text=u"[n]")
        rb1u.place(x=350, y = ys + i*150)
        rb1u.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
        rb2u = tk.Radiobutton(bc_frame,variable=bcvu[i],value=20000+i,text=u"[f]")
        rb2u.place(x=350, y = ys + i*150+30)
        rb2u.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
        rb3u = tk.Radiobutton(bc_frame,variable=bcvu[i],value=30000+i,text=u"[v]")
        rb3u.place(x=350, y = ys + i*150+60)
        rb3u.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
        rb4u = tk.Radiobutton(bc_frame,variable=bcvu[i],value=40000+i,text=u"non")
        rb4u.place(x=350, y = ys + i*150+90)
        rb4u.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        chu1 = tk.Label(bc_frame,text=u"[/m\u00B3]")
        chu1.place(x=530, y = ys + i*150)
        chu1.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
        chu2 = tk.Label(bc_frame,text=u"[/m\u00B2/s]")
        chu2.place(x=530, y = ys + i*150+30)
        chu2.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
        chu3 = tk.Label(bc_frame,text=u"[/m/s]")
        chu3.place(x=530, y = ys + i*150+60)
        chu3.bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        ubc_val1[i] = tk.Entry(bc_frame, width=13)
        if jsp >= 0:
            if Ulabel[jsp] == '1':
                ubc_val1[i].insert(tk.END, Uval[jsp])
        ubc_val1[i].place(x=395,y = ys + i*150 )
        ubc_val1[i].bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        ubc_val2[i] = tk.Entry(bc_frame, width=13)
        if jsp >= 0:
            if Ulabel[jsp] == '2':
                ubc_val2[i].insert(tk.END, Uval[jsp])
        ubc_val2[i].place(x=395,y = ys + i*150+30 )
        ubc_val2[i].bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        ubc_val3[i] = tk.Entry(bc_frame, width=13)
        if jsp >= 0:
            if Ulabel[jsp] == '3':
                ubc_val3[i].insert(tk.END, Uval[jsp])
        ubc_val3[i].place(x=395,y = ys + i*150+60 )
        ubc_val3[i].bind("<MouseWheel>", lambda e:bc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    Done_btn = tk.Button(bc_canvas, font=('',15), text=u'Done')
    Done_btn["command"] = callback_done_boundary_condition_set_window(Planet, dir0, bc_species, bcvl, bcvu, bcset_win, 
                                                                      lbc_val1, lbc_val2, lbc_val3,
                                                                      ubc_val1, ubc_val2, ubc_val3)
    Done_btn.place(x=630, y = 10)

# Calculation setting window
def calculation_set_window(Planet, dir0):

    path = './'+Planet+'/'+dir0+'/settings/calculation_setting.dat'
    dimension_input = ''
    lat_input = ''
    f107_input = ''
    ls_input = ''
    stepmax_input = ''
    sumdt_input = ''
    dtmax_input = ''
    sumdtrot_input = ''
    dtmaxrot_input = ''
    scheme_input = ''
    inversion_input = ''
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            line = f.readlines()
        f.close()

        for i in range(len(line)):
            if 'Dimension' in line[i]:
                dimension_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
            if 'LAT' in line[i]:
                lat_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
            if 'F10.7' in line[i]:
                f107_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
            if 'Ls' in line[i]:
                ls_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
            if 'stepmax' in line[i]:
                stepmax_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
            if 'sumdtst' in line[i]:
                sumdt_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
            if 'dtmaxst' in line[i]:
                dtmax_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
            if 'sumdtrot' in line[i]:
                sumdtrot_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
            if 'dtmaxrot' in line[i]:
                dtmaxrot_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
            if 'scheme' in line[i]:
                scheme_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()
            if 'inversion' in line[i]:
                inversion_input = re.findall(':(.*)',line[i])[0].lstrip().rstrip()

    rbval = {}

    calcset_win = tk.Toplevel()
    calcset_win.title("Calculation settings")
    calcset_win.geometry("720x800")

    calc_canvas = tk.Canvas(calcset_win, width=720,height=800,highlightthickness=0)

    ybar = tk.Scrollbar(calcset_win, orient=tk.VERTICAL) #scroll bar
    ybar.pack(side=tk.RIGHT, fill=tk.Y)
    ybar.config(command=calc_canvas.yview)

    #frame on the main canvas
    calc_frame = tk.Frame(calc_canvas, width=720, height=1400)
    calc_canvas.create_window((0,0), window=calc_frame, anchor=tk.NW, width=calc_canvas.cget('width')) #place frame on the canvas

    calc_canvas.config(yscrollcommand=ybar.set)
    calc_canvas.config(scrollregion=(0,0,720,1400))
    calc_canvas.pack(anchor=tk.NW, expand=1, fill=tk.BOTH)
    calc_canvas.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=0)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '############################################\n')
    text.insert(tk.END, '   Calculation Settings\n')
    text.insert(tk.END, '############################################\n')

    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=100)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '### Dimension setting ###\n')

    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=120)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '   - 1D \n')
    text.insert(tk.END, '   - 2D Lat \n')
    text.insert(tk.END, '   - 2D Rot \n')
    text.insert(tk.END, '   - 3D Rot \n')

    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=100,y=120)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, ': vertical 1D (latitude is fixed. Please set latitude below.)\n')
    text.insert(tk.END, ': latitudenal 2D (LT is fixed to 0 hr)\n')
    text.insert(tk.END, ': Rotational 2D (latitude is fixed. Please set latitude below.)\n')
    text.insert(tk.END, ': Rotational 3D (longitude is fixed. )\n')

    ys = 220

    rbval[1] = tk.IntVar()
    if dimension_input == '':
        dimension_input = 0
    rbval[1].set(int(dimension_input))
    rb1 = tk.Radiobutton(calc_frame,font=("",15),variable=rbval[1],value=1,text=u"1D")
    rb1.place(x=20, y = ys)
    rb1.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    rb2 = tk.Radiobutton(calc_frame,font=("",15),variable=rbval[1],value=2,text=u"2D Lat")
    rb2.place(x=120, y = ys)
    rb2.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    rb3 = tk.Radiobutton(calc_frame,font=("",15),variable=rbval[1],value=3,text=u"2D Rot")
    rb3.place(x=220, y = ys)
    rb3.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    rb4 = tk.Radiobutton(calc_frame,font=("",15),variable=rbval[1],value=4,text=u"3D Rot")
    rb4.place(x=320, y = ys)
    rb4.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=ys+50)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '   - For 1D & 2D Rot mode\n')
    text.insert(tk.END, '      Please set latitude to be calculated. (This is ignored in other modes.)\n')
    text.insert(tk.END, '      Latitude : -90 [deg] ~ 90 [deg] \n')

    LAT = tk.Entry(calc_frame, width=7)
    LAT.insert(tk.END, lat_input)
    LAT.place(x=20,y = ys+120 )
    LAT.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(calc_frame,text=u"[deg]", font=("",15))
    char.place(x=100, y = ys+120)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    ys = 400
    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=ys)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '### Solar EUV F10.7 ###\n')

    char = tk.Label(calc_frame,font=("",15),text=u"F10.7:")
    char.place(x=0, y = ys+50)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    f107 = tk.Entry(calc_frame, width=6)
    f107.insert(tk.END, f107_input)
    f107.place(x=60,y = ys+50 )
    f107.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    ys = 500
    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=ys)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '### Solar longitude [Ls] ###\n')

    char = tk.Label(calc_frame,font=("",15),text=u"Ls:")
    char.place(x=0, y = ys+50)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    ls = tk.Entry(calc_frame, width=6)
    ls.insert(tk.END, ls_input)
    ls.place(x=40,y = ys+50 )
    ls.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(calc_frame,font=("",15),text=u"[deg]")
    char.place(x=120, y = ys+50)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    ys = 600
    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=ys)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '### The number of maximum time steps ###\n')

    char = tk.Label(calc_frame,font=("",15),text=u"Maximum:")
    char.place(x=0, y = ys+50)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    stepmax = tk.Entry(calc_frame, width=15)
    stepmax.insert(tk.END, stepmax_input)
    stepmax.place(x=140,y = ys+50 )
    stepmax.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(calc_frame,font=("",15),text=u"time steps")
    char.place(x=300, y = ys+50)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    ys = 700
    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=ys)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '### Time setting for 1D calculation ###\n')

    char = tk.Label(calc_frame,font=("",15),text=u"Calculate until:")
    char.place(x=0, y = ys+50)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    sumdt = tk.Entry(calc_frame, width=15)
    sumdt.insert(tk.END, sumdt_input)
    sumdt.place(x=140,y = ys+50 )
    sumdt.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(calc_frame,font=("",15),text=u"[sec]")
    char.place(x=300, y = ys+50)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    char = tk.Label(calc_frame,font=("",15),text=u"dt must not excess:")
    char.place(x=0, y = ys+80)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    dtmax = tk.Entry(calc_frame, width=15)
    dtmax.insert(tk.END, dtmax_input)
    dtmax.place(x=140,y = ys+80 )
    dtmax.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(calc_frame,font=("",15),text=u"[sec]")
    char.place(x=300, y = ys+80 )
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    ys = 850
    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=ys)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '### Time setting for 2D & 3D Rotational calculation ###\n')

    char = tk.Label(calc_frame,font=("",15),text=u"Calculate until:")
    char.place(x=0, y = ys+50)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    sumdtrot = tk.Entry(calc_frame, width=15)
    sumdtrot.insert(tk.END, sumdtrot_input)
    sumdtrot.place(x=140,y = ys+50 )
    sumdtrot.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(calc_frame,font=("",15),text=u"[Planetary days]")
    char.place(x=300, y = ys+50)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    char = tk.Label(calc_frame,font=("",15),text=u"dt must not excess:")
    char.place(x=0, y = ys+80)
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    dtmaxrot = tk.Entry(calc_frame, width=15)
    dtmaxrot.insert(tk.END, dtmaxrot_input)
    dtmaxrot.place(x=140,y = ys+80 )
    dtmaxrot.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(calc_frame,font=("",15),text=u"[sec]")
    char.place(x=300, y = ys+80 )
    char.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    ys = 1000
    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=ys)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '### Time advance scheme ###\n')
    text.insert(tk.END, '   - implicit scheme\n')
    text.insert(tk.END, '   - semi-implicit scheme\n')
    text.insert(tk.END, '   - explicit scheme\n')

    rbval[2] = tk.IntVar()
    if scheme_input == '':
        scheme_input = 0
    rbval[2].set(int(scheme_input))
    rb1 = tk.Radiobutton(calc_frame,font=("",15),variable=rbval[2],value=1,text=u"implicit")
    rb1.place(x=20, y = ys+100)
    rb1.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    rb2 = tk.Radiobutton(calc_frame,font=("",15),variable=rbval[2],value=2,text=u"semi-implicit")
    rb2.place(x=120, y = ys+100)
    rb2.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    rb3 = tk.Radiobutton(calc_frame,font=("",15),variable=rbval[2],value=3,text=u"explicit")
    rb3.place(x=260, y = ys+100)
    rb3.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    ys = 1150
    text = tk.Text(calc_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=ys)
    text.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '### Inversion method for implicit scheme ###\n')
    text.insert(tk.END, '   - default\n')

    rbval[3] = tk.IntVar()
    if inversion_input == '':
        inversion_input = 0
    rbval[3].set(int(inversion_input))
    rb1 = tk.Radiobutton(calc_frame,font=("",15),variable=rbval[3],value=1,text=u"default")
    rb1.place(x=20, y = ys+70)
    rb1.bind("<MouseWheel>", lambda e:calc_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))


    # Done button
    Done_btn = tk.Button(calc_canvas, font=('',15), text=u'Done')
    Done_btn["command"] = callback_done_calculation_set_window(Planet, dir0, 
                                                               rbval, 
                                                               LAT, f107, ls, 
                                                               stepmax, 
                                                               sumdt, dtmax, 
                                                               sumdtrot, dtmaxrot, 
                                                               calcset_win)
    Done_btn.place(x=630, y = 10)

# plot window
def plot_window(Planet, dir0):
    plt_win = tk.Toplevel()
    plt_win.title("Plot")
    plt_win.geometry("1000x850")

    plt_canvas = tk.Canvas(plt_win, width=1000,height=850,highlightthickness=0)

    ybar = tk.Scrollbar(plt_win, orient=tk.VERTICAL) #scroll bar
    ybar.pack(side=tk.RIGHT, fill=tk.Y)
    ybar.config(command=plt_canvas.yview)

    species = []

    path = './'+Planet+'/'+dir0+'/settings/plt_species.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            species = f.readlines()

    for isp in range(len(species)):
        species[isp] = species[isp].strip('\n')

    #frame on the main canvas
    ywin = 330 + 30*(len(species) // 5)
    if ywin < 850:
        ywin = 850
    plt_frame = tk.Frame(plt_canvas, width=1000, height=ywin)
    plt_canvas.create_window((0,0), window=plt_frame, anchor=tk.NW, width=plt_canvas.cget('width')) #place frame on the canvas

    plt_canvas.config(yscrollcommand=ybar.set)
    plt_canvas.config(scrollregion=(0,0,1000,ywin))
    plt_canvas.pack(anchor=tk.NW, expand=1, fill=tk.BOTH)
    plt_canvas.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    #plt_frame.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    lines = ['','']
    xrm3 = ['','']
    xrcm3 = ['','']
    xrvmr = ['','']
    yrin = ['','']

    adv = tk.BooleanVar()
    adv.set(True)
    # plot font size
    fs = {}

    path = './'+Planet+'/'+dir0+'/settings/plt_range.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            lines = f.readlines()

    for i in range(len(lines)):
        if 'xrm3' in lines[i]:
            tmp = re.findall('xrm3:(.*)',lines[i])
            xrm3 = tmp[0].split(':')
            xrm3[0].strip('\n').lstrip().rstrip()
            xrm3[1].strip('\n').lstrip().rstrip()
        if 'xrcm3' in lines[i]:
            tmp = re.findall('xrcm3:(.*)',lines[i])
            xrcm3 = tmp[0].split(':')
            xrcm3[0].strip('\n').lstrip().rstrip()
            xrcm3[1].strip('\n').lstrip().rstrip()
        if 'xrvmr' in lines[i]:
            tmp = re.findall('xrvmr:(.*)',lines[i])
            xrvmr = tmp[0].split(':')
            xrvmr[0].strip('\n').lstrip().rstrip()
            xrvmr[1].strip('\n').lstrip().rstrip()
        if 'yr' in lines[i]:
            tmp = re.findall('yr:(.*)',lines[i])
            yrin = tmp[0].split(':')
            yrin[0].strip('\n').lstrip().rstrip()
            yrin[1].strip('\n').lstrip().rstrip()

    #default font size and linewidth of plot
    fs_xlabel = '15'
    fs_ylabel = '15'
    fs_tick   = '15'
    fs_legend = '15'
    linewidth = '1.5'

    path = './'+Planet+'/'+dir0+'/settings/plt_fontsize.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            lines = f.readlines()

    for i in range(len(lines)):
        if 'xlabel' in lines[i]:
            tmp = re.findall('xlabel:(.*)',lines[i])
            fs_xlabel = tmp[0].strip('\n').lstrip().rstrip()
        if 'ylabel' in lines[i]:
            tmp = re.findall('ylabel:(.*)',lines[i])
            fs_ylabel = tmp[0].strip('\n').lstrip().rstrip()
        if 'tick' in lines[i]:
            tmp = re.findall('tick:(.*)',lines[i])
            fs_tick = tmp[0].strip('\n').lstrip().rstrip()
        if 'legend' in lines[i]:
            tmp = re.findall('legend:(.*)',lines[i])
            fs_legend = tmp[0].strip('\n').lstrip().rstrip()
        if 'linewidth' in lines[i]:
            tmp = re.findall('linewidth:(.*)',lines[i])
            linewidth = tmp[0].strip('\n').lstrip().rstrip()

    text = tk.Text(plt_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=0)
    text.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '############################################\n')
    text.insert(tk.END, '   Plot Settings\n')
    text.insert(tk.END, '############################################\n')

    text = tk.Text(plt_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=60)
    text.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '   - Save figure \n')

    text = tk.Text(plt_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=90)
    text.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '         Filename: \n')

    savefname = tk.Entry(plt_frame, width=20)
    savefname.place(x=110,y = 88)
    savefname.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    # save checkbutton
    save_bln = tk.BooleanVar()
    save_bln.set(False)
    save_btn = tk.Checkbutton(plt_frame, anchor="w", variable = save_bln)
    save_btn.place(x = 310, y = 90)
    save_btn.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char_save = tk.Label(plt_frame, width = 30, anchor="w", text = 'Save', font=('', '15'))
    char_save.place(x = 330, y = 90)
    char_save.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    ftrans_bln = tk.BooleanVar()
    ftrans_bln.set(False)
    ftrans_btn = tk.Checkbutton(plt_frame, anchor="w", variable = ftrans_bln)
    ftrans_btn.place(x = 400, y = 90)
    ftrans_btn.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char_ftrans = tk.Label(plt_frame, width = 30, anchor="w", text = 'Background transparent', font=('', '15'))
    char_ftrans.place(x = 420, y = 90)
    char_ftrans.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    text = tk.Text(plt_frame, font=("",15), height=200, width=190, highlightthickness=0)
    text.place(x=0,y=120)
    text.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '   - Select Species \n')

    sp_chk_bln = {}
    for isp in range(len(species)):
        sp_chk_bln[isp] = tk.BooleanVar()
        sp_chk_bln[isp].set(True)
        chk_btn = tk.Checkbutton(plt_frame, anchor="w", variable = sp_chk_bln[isp])
        chk_btn.place(x = 50 + 120*(isp % 5), y = 200 + 30*(isp // 5))
        chk_btn.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

        char_sp = tk.Label(plt_frame, width = 30, anchor="w", text = reaction_unicode(species[isp]), font=('', '15'))
        char_sp.place(x = 75 + 120*(isp % 5), y = 200 + 30*(isp // 5))
        char_sp.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    def all_Select_click():
        for isp in range(len(species)):
            sp_chk_bln[isp].set(True)
    all_Select_Button = tk.Button(plt_frame, text='All Select',command = all_Select_click, font=('', '15'))
    all_Select_Button.place(x=40, y=150)
    all_Select_Button.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    def all_Clear_click():
        for isp in range(len(species)):
            sp_chk_bln[isp].set(False)

    all_Clear_Button = tk.Button(plt_frame, text='All Clear',command = all_Clear_click, font=('', '15'))
    all_Clear_Button.place(x=150, y=150)
    all_Clear_Button.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    def all_Select_ion_click():
        for isp in range(len(species)):
            if '+' in species[isp] or '-' in species[isp]:
                sp_chk_bln[isp].set(True)
    all_Select_ion_Button = tk.Button(plt_frame, text='All Select ion',command = all_Select_ion_click, font=('', '15'))
    all_Select_ion_Button.place(x=250, y=150)
    all_Select_ion_Button.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    def all_Select_neutral_click():
        for isp in range(len(species)):
            if '+' not in species[isp] and '-' not in species[isp]:
                sp_chk_bln[isp].set(True)
    all_Select_neutral_Button = tk.Button(plt_frame, text='All Select neutral',command = all_Select_neutral_click, font=('', '15'))
    all_Select_neutral_Button.place(x=400, y=150)
    all_Select_neutral_Button.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    ys = 30
    rbvar = tk.IntVar()
    rbvar.set(1)

    # Plot mode selection 
    char = tk.Label(plt_frame,text=u"=========== mode selection ===========", font=("",15))
    char.place(x=650, y = ys-30)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    rb1 = tk.Radiobutton(plt_frame,font=("",15),var=rbvar,value=1,text=u"1D stable")
    rb1.place(x=650, y = ys)
    rb1.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    rb2 = tk.Radiobutton(plt_frame,font=("",15),variable=rbvar,value=2,text=u"2D stable")
    rb2.place(x=650, y = ys+30)
    rb2.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"Latitude:", font=("",15))
    char.place(x=650, y = ys+60)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    LATin2 = tk.Entry(plt_frame, width=6)
    LATin2.place(x=750,y = ys+60)
    LATin2.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"(degree)", font=("",15))
    char.place(x=820, y = ys+60)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    rb3 = tk.Radiobutton(plt_frame,font=("",15),variable=rbvar,value=3,text=u"2D & 3D rotation")
    rb3.place(x=650, y = ys+100)
    rb3.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"Latitude:", font=("",15))
    char.place(x=650, y = ys+130)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    LATin3 = tk.Entry(plt_frame, width=6)
    LATin3.place(x=750,y = ys+130)
    LATin3.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"(degree)", font=("",15))
    char.place(x=820, y = ys+130)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"Local Time:", font=("",15))
    char.place(x=650, y = ys+160)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    LTin3 = tk.Entry(plt_frame, width=6)
    LTin3.place(x=750,y = ys+160)
    LTin3.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"(hour)", font=("",15))
    char.place(x=820, y = ys+160)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    # set yrange
    ys = 270
    yr = ['','']
    char = tk.Label(plt_frame,text=u"========== xy range & plot ===========", font=("",15))
    char.place(x=650, y = ys-30)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"yrange:", font=("",15))
    char.place(x=650, y = ys)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    yr[0] = tk.Entry(plt_frame, width=6)
    yr[0].insert(tk.END, yrin[0])
    yr[0].place(x=720,y = ys )
    yr[0].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u" ~ ", font=("",15))
    char.place(x=790, y = ys)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    yr[1] = tk.Entry(plt_frame, width=6)
    yr[1].insert(tk.END, yrin[1])
    yr[1].place(x=815,y = ys )
    yr[1].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u" [km] ", font=("",15))
    char.place(x=885, y = ys)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    # set xrange
    ys = 310
    xr = [['',''],['',''],['','']]
    char = tk.Label(plt_frame,text=u"xrange:", font=("",15))
    char.place(x=650, y = ys)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    xr[0][0]  = tk.Entry(plt_frame, width=6)
    xr[0][0].insert(tk.END, xrm3[0])
    xr[0][0].place(x=720,y = ys )
    xr[0][0].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u" ~ ", font=("",15))
    char.place(x=790, y = ys)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    xr[0][1]  = tk.Entry(plt_frame, width=6)
    xr[0][1].insert(tk.END, xrm3[1])
    xr[0][1].place(x=815,y = ys )
    xr[0][1].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u" [/m\u00B3] ", font=("",15))
    char.place(x=885, y = ys)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    # plot buttons
    ys = 340
    plt_denm_btn = tk.Button(plt_frame, text=u'Plot density [/m\u00B3]', font=('', '15'))
    plt_denm_btn["command"] = callback_plot(Planet, dir0, species, 'density [/m^3]', rbvar, adv, fs, 
                                           sp_chk_bln, '', '', 
                                           xr, yr, LATin2, LATin3, LTin3, 
                                           savefname, save_bln, ftrans_bln)
    plt_denm_btn.place(x=650, y=ys)
    plt_denm_btn.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    char = tk.Label(plt_frame,text=u"xrange:", font=("",15))
    char.place(x=650, y = ys+50)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    xr[1][0]  = tk.Entry(plt_frame, width=6)
    xr[1][0].insert(tk.END, xrcm3[0])
    xr[1][0].place(x=720,y = ys+50 )
    xr[1][0].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u" ~ ", font=("",15))
    char.place(x=790, y = ys+50)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    xr[1][1]  = tk.Entry(plt_frame, width=6)
    xr[1][1].insert(tk.END, xrcm3[1])
    xr[1][1].place(x=815,y = ys+50 )
    xr[1][1].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u" [/cm\u00B3] ", font=("",15))
    char.place(x=885, y = ys+50)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    plt_dencm_btn = tk.Button(plt_frame, text=u'Plot density [/cm\u00B3]', font=('', '15'))
    plt_dencm_btn["command"] = callback_plot(Planet, dir0, species, 'density [/cm^3]', rbvar, adv, fs, 
                                           sp_chk_bln, '', '', 
                                           xr, yr, LATin2, LATin3, LTin3, 
                                           savefname, save_bln, ftrans_bln)
    plt_dencm_btn.place(x=650, y=ys+80)
    plt_dencm_btn.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    ys = 450
    char = tk.Label(plt_frame,text=u"xrange:", font=("",15))
    char.place(x=650, y = ys)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    xr[2][0]  = tk.Entry(plt_frame, width=6)
    xr[2][0].insert(tk.END, xrvmr[0])
    xr[2][0].place(x=720,y = ys )
    xr[2][0].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u" ~ ", font=("",15))
    char.place(x=790, y = ys)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    xr[2][1]  = tk.Entry(plt_frame, width=6)
    xr[2][1].insert(tk.END, xrvmr[1])
    xr[2][1].place(x=815,y = ys )
    xr[2][1].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"  ", font=("",15))
    char.place(x=885, y = ys)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    plt_vmr_btn = tk.Button(plt_frame, text=u'Plot mixing ratio', font=("",15))
    plt_vmr_btn["command"] = callback_plot(Planet, dir0, species, 'mixing ratio', rbvar, adv, fs, 
                                           sp_chk_bln, '', '', 
                                           xr, yr, LATin2, LATin3, LTin3, 
                                           savefname, save_bln, ftrans_bln)
    plt_vmr_btn.place(x=650, y=ys+30)
    plt_vmr_btn.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    # advanced setting
    ys = 530
    char = tk.Label(plt_frame,text=u"========== Advanced setting ==========", font=("",15))
    char.place(x=650, y = ys)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    order_btn = tk.Button(plt_frame, text=u'Colors and orders', font=("",15))
    order_btn["command"] = callback_plot_order(Planet, dir0)
    order_btn.place(x=650, y=ys+30)
    order_btn.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    chk_btn = tk.Checkbutton(plt_frame, anchor="w", variable = adv)
    chk_btn.place(x = 650, y = ys+63)
    chk_btn.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"Apply advanced settings", font=("",15))
    char.place(x=680, y = ys+60)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    # set plot font size
    ys = 630
    char = tk.Label(plt_frame,text=u"============= font size ==============", font=("",15))
    char.place(x=650, y = ys)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"xlabel:", font=("",15))
    char.place(x=650, y = ys+30)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    fs[0] = tk.Entry(plt_frame, width=6)
    fs[0].insert(tk.END,fs_xlabel)
    fs[0].place(x=720,y = ys+30)
    fs[0].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"ylabel:", font=("",15))
    char.place(x=650, y = ys+60)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    fs[1] = tk.Entry(plt_frame, width=6)
    fs[1].insert(tk.END,fs_ylabel)
    fs[1].place(x=720,y = ys+60)
    fs[1].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"ticks:", font=("",15))
    char.place(x=650, y = ys+90)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    fs[2] = tk.Entry(plt_frame, width=6)
    fs[2].insert(tk.END,fs_tick)
    fs[2].place(x=720,y = ys+90)
    fs[2].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"legend:", font=("",15))
    char.place(x=650, y = ys+120)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    fs[3] = tk.Entry(plt_frame, width=6)
    fs[3].insert(tk.END,fs_legend)
    fs[3].place(x=720,y = ys+120)
    fs[3].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u"linewidth:", font=("",15))
    char.place(x=650, y = ys+150)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    fs[4] = tk.Entry(plt_frame, width=6)
    fs[4].insert(tk.END,linewidth)
    fs[4].place(x=720,y = ys+150)
    fs[4].bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    

    # Ratio of 2 species density
    ys = 230 + 30*(len(species) // 5)
    text = tk.Text(plt_frame, font=("",15), height=40, width=60, highlightthickness=0)
    text.place(x=0,y=ys)
    text.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    text.insert(tk.END, '   - Entry 2 Species for density ratio plot\n')

    sp1 = tk.Entry(plt_frame, width=7)
    sp1.place(x=40,y = ys+40 )
    sp1.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    char = tk.Label(plt_frame,text=u" / ", font=("",15))
    char.place(x=120, y = ys+40)
    char.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    sp2 = tk.Entry(plt_frame, width=7)
    sp2.place(x=140,y = ys+40 )
    sp2.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    plt_ratio_btn = tk.Button(plt_frame, text=u'Plot ratio', font=('', '15'))
    plt_ratio_btn["command"] = callback_plot(Planet, dir0, species, 'ratio', rbvar, adv, fs, 
                                             sp_chk_bln, sp1, sp2, 
                                             xr, yr, LATin2, LATin3, LTin3, 
                                             savefname, save_bln, ftrans_bln)
    plt_ratio_btn.place(x=300, y=ys+40)
    plt_ratio_btn.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))  


def plot_order_window(Planet, dir0):
    plt_win = tk.Toplevel()
    plt_win.title("Set plot details")
    plt_win.geometry("1000x800")

    plt_canvas = tk.Canvas(plt_win, width=1000,height=800,highlightthickness=0)

    ybar = tk.Scrollbar(plt_win, orient=tk.VERTICAL) #scroll bar
    ybar.pack(side=tk.RIGHT, fill=tk.Y)
    ybar.config(command=plt_canvas.yview)

    #frame on the main canvas
    lines = []

    path = './'+Planet+'/'+dir0+'/settings/plt_species_order.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            lines = f.readlines()

    detailtext = tk.Text(plt_win, bd = 2, bg = '#EEEEFF', font=('', '15'))
    detailtext.place(x=0,y=100, height=400, width=800 )
    detailtext.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    for i in range(len(lines)):
        if lines[i] != '\n':
            detailtext.insert(tk.END, lines[i])

    plt_detail_btn = tk.Button(plt_win, text=u'Done', font=('', '15'))
    plt_detail_btn["command"] = callback_plot_order_done(Planet, dir0, detailtext, plt_win)
    plt_detail_btn.place(x=400, y=10)
    plt_detail_btn.bind("<MouseWheel>", lambda e:plt_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))


# plot
def plot(Planet, dir0, species, action, adv, fs,
         sp_chk_bln, sp1, sp2, 
         xr, yr, LAT, LT, 
         savelabel, savefname):

    nmax = 1.0e-100
    lines = ['']

    density = [[] for i in range(len(species))]
    mixingratio = [[] for i in range(len(species))]
    altitude = [[] for i in range(len(species))]
    total = []

    fs_xlabel = int(fs[0].get())
    fs_ylabel = int(fs[1].get())
    fs_tick   = int(fs[2].get())
    fs_legend = int(fs[3].get())
    linewidth = float(fs[4].get())

    fname = savefname.get().lstrip().rstrip()

    if savelabel >= 1:
        if fname == '':
            fname = 'Figure.png'

    path = './'+Planet+'/'+dir0+'/settings/plt_species_order.dat'
    if os.path.exists(path) == True:
        with open(path, mode = 'r') as f:
            lines = f.readlines()

    # rule
    # orders are as follows
    # species:Color
    # species{bandle}:Color(#RGB)
    # e.g. 
    #   CO2:#FFFFFF
    #   C3Hn+{C3H+, C3H2+, C3H3+}:#FF9900
    sp_order = ['0' for i in range(len(lines))]
    sp_bandle = ['0' for i in range(len(lines))]
    sp_color = ['0' for i in range(len(lines))]
    for i in range(len(lines)):
        if re.compile('^(?!#)').search(lines[i]):
            if ':' in lines[i]:
                sp_order[i] = re.findall('(.*):',lines[i])[0].lstrip().rstrip()
                sp_color[i] = re.findall(':(.*)',lines[i])[0].lstrip().rstrip()
                if '[' in sp_order[i]:
                    sp_order[i] = re.findall('(.*)\[',lines[i])[0].lstrip().rstrip()
                    sp_bandle[i] = re.findall('\[(.*)\]',lines[i])[0].split(',')
                    for j in range(len(sp_bandle[i])):
                        sp_bandle[i][j] = sp_bandle[i][j].lstrip().rstrip()

    # Applying selected latitude and local time to plot verticacl profiles at a certain location
    if '2D Lat'in action:
        for isp in range(len(species)):
            path1 = './'+Planet+'/'+dir0+'/output/density/2Dstable/'+species[isp]+'.dat'
            if os.path.exists(path1) == True:
                data = np.loadtxt(path1, comments='!')
                ny = len(data[0])-1
                dy = 180.0/(float(ny)-1.0)
                iy = int((float(LAT)+90.0)/dy)
                #print(iy)
                for i in range(len(data)):
                    altitude[isp].append(data[i][0])
                    density[isp].append(data[i][iy+1])
    
    if '3D Rot'in action:
        path1 = './'+Planet+'/'+dir0+'/output/density/3Drot/resolution.dat'
        if os.path.exists(path1) == True:
            data = np.loadtxt(path1, comments='!')
            nx = data[0]
            ny = data[1]
            nz = data[2]
        dx = 24.0/(float(nx)-1.0)
        ix = int((float(LT))/dx)
        if ny >= 2:
            dy = 180.0/(float(ny)-1.0)
            iy = int((float(LAT)+90.0)/dy)
        if ny == 1:
            dy = 1.0
            iy = 0

        for isp in range(len(species)):
            path1 = './'+Planet+'/'+dir0+'/output/density/3Drot/'+str(ix+1)+'/'+species[isp]+'.dat'
            if os.path.exists(path1) == True:
                data = np.loadtxt(path1, comments='!')
                for i in range(len(data)):
                    altitude[isp].append(data[i][0])
                    density[isp].append(data[i][iy+1])

    # Plot density profile [/m^3]
    if 'density [/m^3]' in action:
        if '2D Lat'not in action and '3D Rot'not in action:
            for isp in range(len(species)):
                path = './'+Planet+'/'+dir0+'/output/density/num/'+species[isp]+'.dat'
                if os.path.exists(path) == True:
                    data = np.loadtxt(path, comments='!')
                    for i in range(len(data)):
                        altitude[isp].append(data[i][0])
                        density[isp].append(data[i][1])

        fig = plt.figure(figsize=(8,6))
        ax1 = fig.add_subplot(111)
        x1 = [0]
        y1 = [0]
        label = 0
        if adv.get() == True:
            for i in range(len(sp_order)):
                if sp_order[i] != '0':
                    label = 1
        if label == 0:
            for isp in range(len(species)):
                if species[isp] != 'M' and sp_chk_bln[isp].get() == True:
                    if nmax < np.max(density[isp]):
                        nmax = np.max(density[isp])
                    x1 = density[isp]
                    y1 = altitude[isp]
                    ax1.plot(x1, y1, label=reaction_unicode(species[isp]), linewidth=linewidth)
        if label == 1:
            for i in range(len(sp_order)):
                if sp_bandle[i] == '0':
                    for isp in range(len(species)):
                        if sp_order[i] == species[isp]:
                            if nmax < np.max(density[isp]):
                                nmax = np.max(density[isp])
                            x1 = density[isp]
                            y1 = altitude[isp]
                            ax1.plot(x1, y1, label=reaction_unicode(species[isp]), color = sp_color[i], linewidth=linewidth)
                if sp_bandle[i] != '0':
                    x1 = [0.0 for i in range(len(altitude[0]))]
                    for j in range(len(sp_bandle[i])):
                        for isp in range(len(species)):
                            if sp_bandle[i][j] == species[isp]:
                                if nmax < np.max(density[isp]):
                                    nmax = np.max(density[isp])
                                for k in range(len(x1)):
                                    x1[k] = x1[k] + density[isp][k]
                                y1 = altitude[isp]
                    ax1.plot(x1, y1, label=reaction_unicode(sp_order[i]), color = sp_color[i], linewidth=linewidth)
        ax1.set_xlabel('Density [m'+rf'$^{{{-3}}}$'+']', fontsize=fs_xlabel)
        ax1.set_ylabel('Altitude [km]', fontsize=fs_ylabel)
        plt.tick_params(labelsize=fs_tick)
        plt.rc('legend', fontsize=fs_legend)
        plt.xscale('log')
        xs = xr[0][0].get().lstrip().rstrip()
        xe = xr[0][1].get().lstrip().rstrip()
        if xs == '':
            xs = '1e5'
        if xe == '':
            xe = str(nmax*10)
        xs = float(xs)
        xe = float(xe)
        ax1.set_xlim([xs,xe])

        ys = yr[0].get().lstrip().rstrip()
        ye = yr[1].get().lstrip().rstrip()
        if ys == '':
            ys = np.min(altitude)
        if ye == '':
            ye = np.max(altitude)
        ys = float(ys)
        ye = float(ye)
        ax1.set_ylim([ys,ye])
        plt.legend(loc='best')
        if savelabel == 0:
            plt.show()
        if savelabel == 1:
            plt.savefig('./'+Planet+'/'+dir0+'/figure/'+fname)
        if savelabel == 2:
            plt.savefig('./'+Planet+'/'+dir0+'/figure/'+fname, transparent=True)

    # Plot density profile [/cm^3]
    if 'density [/cm^3]' in action:
        if '2D Lat'not in action and '3D Rot'not in action:
            for isp in range(len(species)):
                path = './'+Planet+'/'+dir0+'/output/density/num/'+species[isp]+'.dat'
                if os.path.exists(path) == True:
                    data = np.loadtxt(path, comments='!')
                    for i in range(len(data)):
                        altitude[isp].append(data[i][0])
                        density[isp].append(data[i][1])

        fig = plt.figure(figsize=(8,6))
        ax2 = fig.add_subplot(111)
        x2 = [0]
        y2 = [0]
        label = 0
        if adv.get() == True:
            for i in range(len(sp_order)):
                if sp_order[i] != '0':
                    label = 1
        if label == 0:
            for isp in range(len(species)):
                if species[isp] != 'M' and sp_chk_bln[isp].get() == True:
                    if nmax < np.max(density[isp]):
                        nmax = np.max(density[isp])
                    x2 = density[isp]
                    for i in range(len(x2)):
                        x2[i] = x2[i] / 1e6
                    y2 = altitude[isp]
                    ax2.plot(x2, y2, label=reaction_unicode(species[isp]), linewidth=linewidth)
        if label == 1:
            for i in range(len(sp_order)):
                if sp_bandle[i] == '0':
                    for isp in range(len(species)):
                        if sp_order[i] == species[isp]:
                            if nmax < np.max(density[isp]):
                                nmax = np.max(density[isp])
                            x2 = density[isp]
                            for j in range(len(x2)):
                                x2[j] = x2[j] / 1e6
                            y2 = altitude[isp]
                            ax2.plot(x2, y2, label=reaction_unicode(species[isp]), color = sp_color[i], linewidth=linewidth)
                if sp_bandle[i] != '0':
                    x2 = [0.0 for i in range(len(altitude[0]))]
                    for j in range(len(sp_bandle[i])):
                        for isp in range(len(species)):
                            if sp_bandle[i][j] == species[isp]:
                                if nmax < np.max(density[isp]):
                                    nmax = np.max(density[isp])
                                for k in range(len(x2)):
                                    x2[k] = x2[k] + density[isp][k]
                                y2 = altitude[isp]
                    for j in range(len(x2)):
                        x2[j] = x2[j] / 1e6
                    ax2.plot(x2, y2, label=reaction_unicode(sp_order[i]), color = sp_color[i], linewidth=linewidth)

        ax2.set_xlabel('Density [cm'+rf'$^{{{-3}}}$'+']', fontsize=fs_xlabel)
        ax2.set_ylabel('Altitude [km]', fontsize=fs_ylabel)
        plt.tick_params(labelsize=fs_tick)
        plt.rc('legend', fontsize=fs_legend)
        plt.xscale('log')
        xs = xr[1][0].get()
        xe = xr[1][1].get()
        if xs == '':
            xs = '1e-1'
        if xe == '':
            xe = str(nmax*10/1e6)
        xs = float(xs)
        xe = float(xe)
        ax2.set_xlim([xs,xe])

        ys = yr[0].get()
        ye = yr[1].get()
        if ys == '':
            ys = np.min(altitude)
        if ye == '':
            ye = np.max(altitude)
        ys = float(ys)
        ye = float(ye)
        ax2.set_ylim([ys,ye])
        plt.legend(loc='best')
        if savelabel == 0:
            plt.show()
        if savelabel == 1:
            plt.savefig('./'+Planet+'/'+dir0+'/figure/'+fname)
        if savelabel == 2:
            plt.savefig('./'+Planet+'/'+dir0+'/figure/'+fname, transparent=True)

    # Plot mixing ratio profile [mol/mol]
    if 'mixing ratio' in action:
        if '2D Lat'not in action and '3D Rot'not in action:
            for isp in range(len(species)):
                path = './'+Planet+'/'+dir0+'/output/density/vmr/vmr_'+species[isp]+'.dat'
                if os.path.exists(path) == True:
                    data = np.loadtxt(path, comments='!')
                    for i in range(len(data)):
                        altitude[isp].append(data[i][0])
                        mixingratio[isp].append(data[i][1])
        
        if '2D Lat' in action or '3D Rot' in action:
            for iz in range(len(altitude[0])):
                tmp = 0
                for isp in range(len(species)):
                    if species[isp] != 'M':
                        tmp = tmp + density[isp][iz]
                total.append(tmp)
        
        for isp in range(len(species)):
            for iz in range(len(density[isp])):
                mixingratio[isp].append(density[isp][iz]/total[iz])

        fig = plt.figure(figsize=(8,6))
        ax3 = fig.add_subplot(111)
        x3 = [0]
        y3 = [0]
        for isp in range(len(species)):
            if species[isp] != 'M' and sp_chk_bln[isp].get() == True:
                x3 = mixingratio[isp]
                y3 = altitude[isp]
                ax3.plot(x3, y3, label=species[isp], linewidth=linewidth)
        ax3.set_xlabel('Volume mixing ratio', fontsize=fs_xlabel)
        ax3.set_ylabel('Altitude [km]', fontsize=fs_ylabel)
        plt.tick_params(labelsize=fs_tick)
        plt.rc('legend', fontsize=fs_legend)
        plt.xscale('log')
        xs = xr[2][0].get().lstrip().rstrip()
        xe = xr[2][1].get().lstrip().rstrip()
        if xs == '':
            xs = '1e-15'
        if xe == '':
            xe = '1'
        xs = float(xs)
        xe = float(xe)
        ax3.set_xlim([xs,xe])

        ys = yr[0].get()
        ye = yr[1].get()
        if ys == '':
            ys = np.min(altitude)
        if ye == '':
            ye = np.max(altitude)
        ys = float(ys)
        ye = float(ye)
        ax3.set_ylim([ys,ye])
        plt.legend(loc='upper left')
        if savelabel == 0:
            plt.show()
        if savelabel == 1:
            plt.savefig('./'+Planet+'/'+dir0+'/figure/'+fname)
        if savelabel == 2:
            plt.savefig('./'+Planet+'/'+dir0+'/figure/'+fname, transparent=True)

    # Plot 2 species  ratio profile [mol/mol]
    if 'ratio' in action and 'mixing' not in action:
        if '2D Lat'not in action and '3D Rot'not in action:
            for isp in range(len(species)):
                path = './'+Planet+'/'+dir0+'/output/density/num/'+species[isp]+'.dat'
                if os.path.exists(path) == True:
                    data = np.loadtxt(path, comments='!')
                    for i in range(len(data)):
                        altitude[isp].append(data[i][0])
                        density[isp].append(data[i][1])

        fig = plt.figure(figsize=(8,6))
        ax4 = fig.add_subplot(111)
        x4 = [0]
        y4 = [0]
        for isp in range(len(species)):
            csp2 = sp2.get().lstrip().rstrip()
            if species[isp] == csp2:
                x4 = density[isp]
                y4 = altitude[isp]
            unicsp2 = reaction_unicode(csp2)
        for isp in range(len(species)):
            csp1 = sp1.get().lstrip().rstrip()
            if species[isp] == csp1:
                x5 = density[isp]
                for i in range(len(x5)):
                    x5[i] = x5[i] / x4[i]
                y5 = altitude[isp]
                unicsp1 = reaction_unicode(csp1)
                ax4.plot(x5, y4, label=unicsp1+' / '+unicsp2+' ratio', linewidth=linewidth)
        ax4.set_xlabel(unicsp1+' / '+unicsp2+' ratio', fontsize=fs_xlabel)
        ax4.set_ylabel('Altitude [km]', fontsize=fs_ylabel)
        plt.tick_params(labelsize=fs_tick)
        plt.rc('legend', fontsize=fs_legend)
        ys = yr[0].get()
        ye = yr[1].get()
        if ys == '':
            ys = min(altitude)
        if ye == '':
            ye = max(altitude)
        ys = float(ys)
        ye = float(ye)
        ax4.set_ylim([ys,ye])
        plt.xscale('log')
        plt.legend(loc='upper left')
        if savelabel == 0:
            plt.show()
        if savelabel == 1:
            plt.savefig('./'+Planet+'/'+dir0+'/figure/'+fname)
        if savelabel == 2:
            plt.savefig('./'+Planet+'/'+dir0+'/figure/'+fname, transparent=True)


# detailed reference window if no doi link is available
def ref_window(ref, ref_info):
    ref_win = tk.Toplevel() #display on the main window
    ref_win.title("Detailed Reference Information")
    ref_win.geometry("600x800")

    text = tk.Text(ref_win, font=("",15), height=200, width=190, highlightthickness=0)
    text.pack()

    text.insert(tk.END, '############################################\n')
    text.insert(tk.END, ' Detailed Reference Information\n')
    text.insert(tk.END, '############################################\n')
    text.insert(tk.END, '\n')
    text.insert(tk.END, ' Reference:   '+ref+'\n')
    text.insert(tk.END, '\n')
    text.insert(tk.END, ' Detail:   '+ref_info+'\n')

# insert hyperlink tkinter windows
class tkHyperlink:

    def __init__(self, text):

        self.text = text
        self.text.tag_config("hyper", foreground="blue", underline=1)
        self.text.tag_bind("hyper", "<Button-1>", self._click)

        self.reset()

    def reset(self):
        self.links = {}
        self.url = {}

    def add(self, action, url):
        tag = "hyper-%d" % len(self.links)
        self.links[tag] = action
        self.url[tag] = url
        return "hyper", tag

    def _click(self, event):
        for tag in self.text.tag_names(CURRENT):
            if tag[:6] == "hyper-":
                self.links[tag](self.url[tag])
                return

# Make Reaction & rate list window
def reaction_window(iplnt, Planet, list_s, list_e, dir0, version,
                    reaction_chk_bln, fix_species_bln, input_species_char,
                    search_list):

    #open a web browser when a hyperlink on a reference is clicked
    def callback_web(ref):
        for i in range(len(doi_list)):
            doi_list[i][0] = doi_list[i][0].lstrip().rstrip()
            if doi_list[i][0] in ref:
                ref_window(doi_list[i][0], doi_list[i][1])
                if 'doi:' in doi_list[i][1]:
                    doi = re.findall('doi:(.*)', doi_list[i][1])
                    doi = doi[0].lstrip().rstrip().rstrip('.')
                    url = 'https://doi.org/'+doi
                    webbrowser.open(url)

    #read saved reaction list
    global reaction_win_label
    if reaction_win_label == 0:
        loaded_reaction_list=[]
        path = './'+Planet+'/'+dir0+'/settings/reaction_list.dat'
        if os.path.exists(path) == True:
            with open(path, mode = 'r') as f:
                loaded_reaction_list = f.readlines()
            f.close()

        for jch in range(len(loaded_reaction_list)):
            # Treatment for old expressions for Tn, Ti, Te
            if 'T(neutral)' in loaded_reaction_list[jch]:
                loaded_reaction_list[jch] = re.sub('T\(neutral\)','Tn',loaded_reaction_list[jch])
            if 'T(ion)' in loaded_reaction_list[jch]:
                loaded_reaction_list[jch] = re.sub('T\(ion\)','Ti',loaded_reaction_list[jch])
            if 'T(electron)' in loaded_reaction_list[jch]:
                loaded_reaction_list[jch] = re.sub('T\(electron\)','Te',loaded_reaction_list[jch])
            label=0
            for ich in range(len(reaction_rate_list)):
                if ich >= list_s and ich <= list_e: #reaction list [start ~ end]
                    if loaded_reaction_list[jch].strip('\n') == reaction_rate_list[ich]:
                        label=1
                        reaction_chk_bln[ich].set(True)
            if label==0:
                print('following reaction is not seen in the reaction list.\n\n')
                print(loaded_reaction_list[jch])
        reaction_win_label += 1

    # frame, scrollbar
    nlist = 0 #height of the window frame
    rate_eq = [] #

    #make a list of rate coefficient
    for ich in range(len(reaction_rate_list)):
        if ich >= list_s and ich <= list_e: #reaction list [start ~ end]

            if search_list[ich] == 1:
                nlist = nlist + 1

                rate_eq_tmp = re.findall(':(.*)', reaction_rate_list[ich]) #take reaction coefficient and/or reference
                nlist = nlist + len(rate_eq_tmp[0].split('&&')) - 1 #if this coefficient has 2 cases. ex) Temperature dependence

    ### set each canvas ########################################################
    # scroll bar
    ybar = tk.Scrollbar(win, orient=tk.VERTICAL) #scroll bar
    xbar = tk.Scrollbar(win, orient=tk.HORIZONTAL) #scroll bar
    xbar.pack(side=tk.BOTTOM, fill=tk.X)
    ybar.pack(side=tk.RIGHT, fill=tk.Y)
    #upper canvas
    upper_canvas = tk.Canvas(win, width=100,height=130, highlightthickness=0)
    upper_canvas.pack(side=tk.TOP, fill=tk.BOTH)
    #main canvas
    main_canvas = tk.Canvas(win, width=1280,height=600, highlightthickness=0)
    main_canvas.pack(anchor=tk.NW, expand=1, fill=tk.BOTH)
    main_canvas.bind('<MouseWheel>',lambda e:main_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    #frame on the main canvas
    frame = tk.Frame(main_canvas, width=1280,height=nlist*30+120)
    main_canvas.create_window((0,0), window=frame, anchor=tk.NW, width=main_canvas.cget('width')) #place frame on the canvas
    frame.bind('<MouseWheel>',lambda e:main_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    #lower canvas
    lower_canvas = tk.Canvas(win, width=1280,height=60, highlightthickness=0)
    lower_canvas.pack(side=tk.BOTTOM, fill=tk.BOTH)

    txt = Planet + '    |  Curent project:   ./'+Planet+'/'+dir0
    labelTop01 = tk.Label(upper_canvas, font=('',25), text=txt)
    labelTop01.place(x=20,y=0)
    txt = 'Version:  '+str(version)
    labelTop02 = tk.Label(upper_canvas, font=('',25), text=txt)
    labelTop02.place(x=620,y=0)
    labelTop1 = tk.Label(upper_canvas, text=u'                        Chemical Reactions')
    labelTop1.place(x=0,y=90)
    labelTop2 = tk.Label(upper_canvas, text=u'|      Rates')
    labelTop2.place(x=400, y=90)
    labelTop3 = tk.Label(upper_canvas, text=u'|              Label ')
    labelTop3.place(x=780, y=90)
    labelTop4 = tk.Label(upper_canvas, text=u'|         Reference       ')
    labelTop4.place(x=1080, y=90)

    hbar = tk.Label(upper_canvas, text=u'------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
    hbar.place(x=0,y=110)

    #set scroll bar
    #config function can access attributes of object
    xbar.config(command=main_canvas.xview)
    ybar.config(command=main_canvas.yview)
    main_canvas.config(xscrollcommand=xbar.set, yscrollcommand=ybar.set)
    main_canvas.config(scrollregion=(0,0,1280,nlist*30+120)) #(upper_left_x,upper_left_y,lower_right_x,lower_right_y)

    ########## Buttons ############

    ### search ###
    #search input text
    search_text = tk.Entry(upper_canvas, width=20)
    search_text.place(x=20,y=36)
    # set search button
    Search_btn = tk.Button(upper_canvas, text=u'Search reaction')
    Search_btn["command"] = callback_search_reaction_list(main_canvas, xbar, ybar, frame, upper_canvas, lower_canvas,
                                                          iplnt, Planet, list_s, list_e, dir0, version,
                                                          reaction_chk_bln, fix_species_bln, input_species_char,
                                                          search_species, search_text, search_list)
    Search_btn.place(x=220, y=40)

    ### update button ###
    # TO DO "update" name should be changed later
    Update_btn = tk.Button(upper_canvas, text=u'Update list')
    Update_btn["command"] = callback_update_reaction_window(main_canvas, xbar, ybar, frame, upper_canvas, lower_canvas,
                                                            iplnt, Planet, list_s, list_e, dir0, version,
                                                            reaction_chk_bln, fix_species_bln, input_species_char,
                                                            search_list)
    Update_btn.place(x=220, y=65)

    ### grid button ###
    Grid_btn = tk.Button(upper_canvas, text=u'Set z Grid')
    Grid_btn["command"] = callback_grid_window(Planet, dir0)
    Grid_btn.place(x=500, y=40)

    ### input and fix button ###
    Fix_btn = tk.Button(upper_canvas, text=u'Set input species')
    Fix_btn["command"] = callback_input_window(Planet, dir0, fix_species_bln, input_species_char)
    Fix_btn.place(x=350, y=40)

    ### boundary condition button ###
    bc_btn = tk.Button(upper_canvas, text=u'Boundary Condition')
    bc_btn["command"] = callback_boundary_condition_input_window(Planet, dir0)
    bc_btn.place(x=350, y=65)

    ### calculation setting button ###
    calc_btn = tk.Button(upper_canvas, text=u'Calculation settings')
    calc_btn["command"] = callback_calculation_set_window(Planet, dir0)
    calc_btn.place(x=500, y=65)

    ### help button ###
    help_btn = tk.Button(upper_canvas, text=u'Help')
    help_btn["command"] = callback_help_window()
    help_btn.place(x=680, y=40)

    ### done button ###
    exit_btn = tk.Button(upper_canvas, text=u'Exit')
    exit_btn["command"] = callback_exit_window(win)
    exit_btn.place(x=750, y=40)
    ############################################################################

    ### main canvas ############################################################
    # display Chemical reactions and Check Buttons
    ilist = 0
    all_merged_reaction_L_text = ""
    all_merged_array_text = ""
    all_merged_reaction_R_text = ""
    all_merged_rate_text = ""

    #label
    chk_label_all = tk.Text(frame, font=("",12), width = 50, height=nlist*30, highlightthickness=0)
    chk_label_all.place(x = 810, y = 12)
    hyperlink_label = tkHyperlink(chk_label_all)
    chk_label_all.bind("<MouseWheel>", lambda e:main_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    #reference
    chk_reference_all = tk.Text(frame, font=("",12), width = 100, height=nlist*30, highlightthickness=0)
    chk_reference_all.place(x = 1100, y = 12)
    hyperlink_reference = tkHyperlink(chk_reference_all)
    chk_reference_all.bind("<MouseWheel>", lambda e:main_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))

    reaction_chk_bln_old = [-1 for i in range(len(reaction_rate_list))]
    for ich in range(len(reaction_rate_list)):
        if ich >= list_s and ich <= list_e and search_list[ich] == 1:
            if reaction_chk_bln[ich].get()==True:
                reaction_chk_bln_old[ich] = 1
            if reaction_chk_bln[ich].get()==False:
                reaction_chk_bln_old[ich] = 0

    for ich in range(len(reaction_rate_list)):
        if ich >= list_s and ich <= list_e and search_list[ich] == 1:
            #print(reaction_rate_list[ich])
            # replace ^N, _N, ^+, ^- >> unicode

            # reaction display
            reaction_find = re.findall('(.*):', reaction_rate_list[ich]) #species part
            reaction = reaction_find[0]
            reaction = reaction_unicode(reaction)
            reaction = re.sub('\s+', '  ', reaction)
            
            #divide reaction species into reactants(left side) and products(right side)
            reaction_Ldisplay = re.findall('(.*)\u2192', reaction) #\u2192 = arrow
            reaction_Rdisplay = re.findall('\u2192(.*)', reaction)
            #remove spaces on the left and right
            reaction_Ldisplay[0] = reaction_Ldisplay[0].lstrip().rstrip()
            reaction_Rdisplay[0] = reaction_Rdisplay[0].lstrip().rstrip()

            # @ denotes reference. Some reaction have reference with coefficient
            rate_ref_label_tmp = re.findall(':(.*)', reaction_rate_list[ich])
            rate_ref_label_all_cases = rate_ref_label_tmp[0].split('&&')
            rate_all_cases = ['' for i in range(len(rate_ref_label_all_cases))]
            ref_all_cases = ['' for i in range(len(rate_ref_label_all_cases))]
            label_all_cases = ['' for i in range(len(rate_ref_label_all_cases))]
            for i in range(len(rate_ref_label_all_cases)):
                if '@' not in rate_ref_label_all_cases[i] and '#' not in rate_ref_label_all_cases[i]:
                    rate_all_cases[i] = rate_ref_label_all_cases[i]
                if '@' in rate_ref_label_all_cases[i] and '#' not in rate_ref_label_all_cases[i]:
                    rate_tmp = re.findall('(.*)@', rate_ref_label_all_cases[i])
                    rate_all_cases[i] = rate_tmp[0]
                    ref_tmp  = re.findall('@(.*)', rate_ref_label_all_cases[i])
                    ref_all_cases[i] = ref_tmp[0]
                if '@' not in rate_ref_label_all_cases[i] and '#' in rate_ref_label_all_cases[i]:
                    rate_tmp = re.findall('(.*)#', rate_ref_label_all_cases[i])
                    rate_all_cases[i] = rate_tmp[0]
                    label_tmp  = re.findall('#(.*)', rate_ref_label_all_cases[i])
                    label_all_cases[i] = label_tmp[0]
                if '@' in rate_ref_label_all_cases[i] and '#' in rate_ref_label_all_cases[i]:
                    rate_tmp = re.findall('(.*)@', rate_ref_label_all_cases[i])
                    if '#' in rate_tmp[0]:
                        rate_tmp1 = re.findall('(.*)#', rate_tmp[0])
                        rate_all_cases[i] = rate_tmp1[0]
                        ref_tmp  = re.findall('@(.*)', rate_ref_label_all_cases[i])
                        ref_all_cases[i] = ref_tmp[0]
                        label_tmp  = re.findall('#(.*)@', rate_ref_label_all_cases[i])
                        label_all_cases[i] = label_tmp[0]
                    if '#' not in rate_tmp[0]:
                        rate_all_cases[i] = rate_tmp[0]
                        ref_tmp  = re.findall('@(.*)#', rate_ref_label_all_cases[i])
                        ref_all_cases[i] = ref_tmp[0]
                        label_tmp  = re.findall('#(.*)', rate_ref_label_all_cases[i])
                        label_all_cases[i] = label_tmp[0]
                rate_all_cases[i] = rate_unicode(rate_all_cases[i])
                rate_all_cases[i] = rate_all_cases[i].lstrip().rstrip() #remove spaces on the left and right
                ref_all_cases[i] = ref_all_cases[i].lstrip().rstrip()
                label_all_cases[i] = label_all_cases[i].lstrip().rstrip()

            ### check buttons
            #anchor w means west side of the center
            chk_reaction_btn = tk.Checkbutton(frame, width = 30, anchor="w", variable = reaction_chk_bln[ich])
            chk_reaction_btn.place(x = 5, y = 10 + (ilist * 30))

            ### reactions
            all_merged_reaction_L_text += reaction_Ldisplay[0]
            all_merged_array_text += "\u2192"
            all_merged_reaction_R_text += reaction_Rdisplay[0]

            ### rates
            #some reactions have multiple coefficients, not only 1!!
            for i in range(len(rate_all_cases)):
                ilist = ilist + 1
                all_merged_reaction_L_text += "\n\n"
                all_merged_array_text += "\n\n"
                all_merged_reaction_R_text += "\n\n"
                all_merged_rate_text += rate_all_cases[i] + "\n\n"
                chk_reference_all.insert(tk.INSERT, ref_all_cases[i], hyperlink_reference.add(callback_web, ref_all_cases[i]))
                chk_reference_all.insert(tk.INSERT, "\n\n")
                chk_label_all.insert(tk.INSERT, label_all_cases[i], hyperlink_label.add(callback_web, label_all_cases[i]))
                chk_label_all.insert(tk.INSERT, "\n\n")
            ilist = ilist - 1


            ilist = ilist + 1
    ############################################################################

    ### place Text object on the frame
    #reaction
    chk_reaction_L = tk.Text(frame, font=("",12), width = 45, height=nlist*30, highlightthickness=0)
    chk_reaction_L.insert(tk.INSERT, all_merged_reaction_L_text)
    chk_reaction_L.place(x = 25, y = 12)
    chk_reaction_L.bind("<MouseWheel>", lambda e:main_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    chk_array = tk.Text(frame, font=("",12), width = 45, height=nlist*30, highlightthickness=0)
    chk_array.insert(tk.INSERT, all_merged_array_text)
    chk_array.place(x = 150, y = 12)
    chk_array.bind("<MouseWheel>", lambda e:main_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    chk_reaction_R = tk.Text(frame, font=("",12), width = 45, height=nlist*30, highlightthickness=0)
    chk_reaction_R.insert(tk.INSERT, all_merged_reaction_R_text)
    chk_reaction_R.place(x = 200, y = 12)
    chk_reaction_R.bind("<MouseWheel>", lambda e:main_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))
    #rate
    chk_rate_all = tk.Text(frame, font=("",12), width = 50, height=nlist*30, highlightthickness=0)
    chk_rate_all.insert(tk.INSERT, all_merged_rate_text)
    chk_rate_all.place(x = 400, y = 12)
    chk_rate_all.bind("<MouseWheel>", lambda e:main_canvas.yview_scroll(-1*(1 if e.delta>0 else -1),'units'))


    ### Lower canvas ###########################################################
    # Buttons and labels
    bbar = tk.Label(lower_canvas, text=u'------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
    bbar.place(x=0,y=0)

    Analysis_btn = tk.Button(lower_canvas, text=u'Reaction Analysis')
    Analysis_btn["command"] = callback_reaction_analysis('', iplnt, reaction_chk_bln, fix_species_bln, dir0)
    Analysis_btn.place(x=350, y=20)

    Output_btn = tk.Button(lower_canvas, text=u'Output f90 module')
    Output_btn["command"] = callback_reaction_analysis('Output', iplnt, reaction_chk_bln, fix_species_bln, dir0)
    Output_btn.place(x=500, y=20)

    run_btn = tk.Button(lower_canvas, text=u'Output f90 module\n& Run model')
    run_btn["command"] = callback_reaction_analysis('Run', iplnt, reaction_chk_bln, fix_species_bln, dir0)
    run_btn.place(x=650, y=20)

    plot_btn = tk.Button(upper_canvas, text=u'Plot setting', font=('', '15'))
    plot_btn["command"] = callback_plot_window(Planet, dir0)
    plot_btn.place(x=680, y=65)

    # set all the check buttons in displayed window as True
    def all_Select_click():
        for ich in range(len(reaction_rate_list)):
            if ich >= list_s and ich <= list_e and search_list[ich] == 1:
                reaction_chk_bln[ich].set(True)
    all_Select_Button = tk.Button(lower_canvas, text='All Select',command = all_Select_click)
    all_Select_Button.place(x=10, y=20)

    # set all the check buttons in displayed window as False
    def all_Clear_click():
        for ich in range(len(reaction_rate_list)):
            if ich >= list_s and ich <= list_e and search_list[ich] == 1:
                reaction_chk_bln[ich].set(False)
    all_Clear_Button = tk.Button(lower_canvas, text='All Clear',command = all_Clear_click)
    all_Clear_Button.place(x=80, y=20)

    # undo check buttons in displayed window 
    reaction_chk_bln_new = [-1 for i in range(len(reaction_rate_list))]
    def undo_click():
        for ich in range(len(reaction_rate_list)):
            if ich >= list_s and ich <= list_e and search_list[ich] == 1:
                if reaction_chk_bln[ich].get()==True:
                    reaction_chk_bln_new[ich] = 1
                if reaction_chk_bln[ich].get()==False:
                    reaction_chk_bln_new[ich] = 0
        for ich in range(len(reaction_rate_list)):
            if ich >= list_s and ich <= list_e and search_list[ich] == 1:
                if reaction_chk_bln_old[ich] == 1:
                    reaction_chk_bln[ich].set(True)
                if reaction_chk_bln_old[ich] == 0:
                    reaction_chk_bln[ich].set(False)
    undo_Button = tk.Button(lower_canvas, text='Undo',command = undo_click)
    undo_Button.place(x=150, y=20)

    # redo check buttons in displayed window: it works after undo button pressed, 
    # but does nothing without pressing undo button before
    def redo_click():
        for ich in range(len(reaction_rate_list)):
            if ich >= list_s and ich <= list_e and search_list[ich] == 1:
                if reaction_chk_bln_new[ich] == 1:
                    reaction_chk_bln[ich].set(True)
                if reaction_chk_bln_new[ich] == 0:
                    reaction_chk_bln[ich].set(False)
    redo_Button = tk.Button(lower_canvas, text='Redo',command = redo_click)
    redo_Button.place(x=200, y=20)
    ############################################################################

##### Planet Select Window #####################################################
#header
Select_planets = tk.Label(win, font=('', '20'), text=u'Choose Planet')
#Select_planets.pack()
Select_planets.place(x=30, y=10) #place label object
hbar = tk.Label(win, text=u'----------------------------------------------------')
#hbar.pack()
hbar.place(x=0,y=50)

#planet button
Planet_btn = {}
for iplnt in range(len(Planet_list)):
    if Planet_list[iplnt][0] != 'END':
        Planet_btn[iplnt] = tk.Button(win, font=('', '20'), text=Planet_list[iplnt][0]) #text=planet name
        Planet_btn[iplnt]["command"] = callback_directory_window(iplnt, list_s, list_e,
                                                                reaction_chk_bln, fix_species_bln, input_species_char,
                                                                search_list) #display reaction window if pushed
        #Planet_btn[iplnt].pack()
        Planet_btn[iplnt].place(x=30, y=90+iplnt*60)

#help button
help_btn = tk.Button(win, font=('', '20'), text=u'Help')
help_btn["command"] = callback_help_window()
#help_btn.pack()
help_btn.place(x=850, y=10)
done_btn = tk.Button(win, font=('', '20'), text=u'Exit')
done_btn["command"] = callback_exit_window(win)
#done_btn.pack()
done_btn.place(x=1050, y=10)
################################################################################

# run main window
win.mainloop()
