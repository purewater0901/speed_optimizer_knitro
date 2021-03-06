import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('../position_result.csv')

data.plot(x='time', figsize=(20,18))

plt.ylim([-2.0, 5.0])
plt.xlabel('time[s]')
plt.ylabel('speed[m/s]')
plt.title('time optimized speed')
plt.savefig('../result/time_optimized_speed.png')
plt.show()