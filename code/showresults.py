import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

sample_info = pd.read_csv("sample_info.csv")

plt.plot(sample_info["scalars[0]"], sample_info["scalars[1]"],
         "b.", markersize=1, alpha=0.2)
plt.show()

