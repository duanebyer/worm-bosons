#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import numpy as np

""" Project 1, part 2 """
""" main.cpp compiled to /main """

""" Parameters """ 
num_events = 100000
num_burn = 10000
N = 2
mu = 1.4
beta = 12
t = 1
epsilon_list = [0.03,0.02,0.01,0.008,0.005,0.002,0.001]
energy_list = np.array([])
number_list = np.array([])

""" Running worm code and saving results"""
for epsilon in epsilon_list:
    process = subprocess.run(['./main', str(num_events), str(num_burn), str(N), str(epsilon), str(mu), str(beta), str(t)],  
               capture_output=True, text=True)
    lines = process.stdout.splitlines()
    energy_list = np.append(energy_list, float(lines[-2][12:]))
    number_list = np.append(number_list, float(lines[-1][12:]))
    
""" Polynomial fit """
energy_fit = np.polyfit(epsilon_list, energy_list, deg=3)
number_fit = np.polyfit(epsilon_list, number_list, deg=3)

print("<n> =", number_fit[-1])
print("<e> =", energy_fit[-1])