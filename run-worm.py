#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import multiprocessing as mp
import numpy as np

num_events = 2**18
num_burn = 10000
N_axis = np.array([20])
epsilon_axis = np.array([0.01])
mu_axis = np.array([1.0])
beta_axis = np.linspace(0.75, 1.1, 40)
t_axis = np.array([1.])

N, epsilon, mu, beta, t = np.meshgrid(N_axis, epsilon_axis, mu_axis, beta_axis, t_axis, indexing='ij')
N = np.ravel(N)
epsilon = np.ravel(epsilon)
mu = np.ravel(mu)
beta = np.ravel(beta)
t = np.ravel(t)

def process(params):
    N, epsilon, mu, beta, t = params
    process = subprocess.run(['./worm', str(num_events), str(num_burn), str(int(N)), str(epsilon), str(mu), str(beta), str(t)], capture_output=True, text=True)
    lines = process.stdout.splitlines()
    if process.returncode != 0:
        raise "Worm algorithm had an error"
    result = []
    for line in lines[-3:]:
        tokens = line.split(' ')
        result.append([float(tokens[1]), float(tokens[3])])
    return result

with mp.Pool(8) as p:
    params = np.column_stack((N, epsilon, mu, beta, t))
    outputs = p.map(process, params)
    outputs = np.reshape(np.array(outputs), (np.size(N), 6))
    outputs = np.column_stack((params, outputs))
    np.savetxt('worm-output.csv', outputs, header='N ε μ β t <e> <e>_err <n> <n>_err χ χ_err', comments='')

