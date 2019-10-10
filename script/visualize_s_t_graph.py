import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('../position_result.csv')

data.plot(x='position', y='time')
for i in range(-10,10):
    plt.vlines(20+i*0.01, 0, 20, colors='black')
    plt.vlines(20+i*0.01, 80, 90, colors='black')

plt.xlabel('position[m]')
plt.ylabel('time[s]')
plt.title('time optimized speed')
plt.savefig('../result/st_graph.png')
plt.show()