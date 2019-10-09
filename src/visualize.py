import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('../position_result.csv')

data.plot(x='position')

plt.ylim([0, 10])
plt.xlabel('position[m]')
plt.ylabel('speed[m/s]')
plt.title('time optimized speed')
plt.savefig('../result/time_optimized_speed.png')
plt.show()