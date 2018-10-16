from matplotlib import pyplot as plt
import numpy as np

datatype = np.dtype([
    ('N', int), 
    ('y0', float), 
    ('P', float), 
    ('tau0', float), 
    ('dt', float),
    ('umax', float)])
data = np.loadtxt('bingham_poiseuille.dat', datatype, skiprows=1)
print(data)

num_N_vals = len(np.unique(data['N']))
# Set up figure environment
fig = plt.figure(num=num_N_vals, figsize=(6,6))
fig.clf()

for i, N in enumerate(np.unique(data['N'])):
    ax = fig.add_subplot(num_N_vals,1,i+1)
    mask = data['N'] == N
    for datapoint in data[mask]:
        if datapoint['umax'] == 0.0:
            color = 'red'
        else:
            color = 'green'
        plt.plot(datapoint['tau0'],  datapoint['dt'], marker='.', markerfacecolor=None, color=color)

plt.show()

num_tau_vals = len(np.unique(data['tau0']))
# Set up figure environment
fig = plt.figure(num=num_tau_vals, figsize=(3,10))
fig.clf()

for i, tau in enumerate(np.unique(data['tau0'])):
    ax = fig.add_subplot(num_tau_vals,1,i+1)
    mask = data['tau0'] == tau
    for datapoint in data[mask]:
        if datapoint['umax'] == 0.0:
            color = 'red'
        else:
            color = 'green'
        plt.plot(datapoint['N'],  datapoint['dt'], marker='.', markerfacecolor=None, color=color)

plt.show()
