import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Read CSV using pandas
df0 = pd.read_csv("../result/data/ranker_socket_stream_coop1_0.csv")
df1 = pd.read_csv("../result/data/ranker_socket_stream_coop1_1.csv")
df2 = pd.read_csv("../result/data/ranker_socket_stream_coop1_2.csv")
df3 = pd.read_csv("../result/data/ranker_socket_stream_coop1_3.csv")
    
# Extract values from dataframe
x0 = df0['AccessFrequency(line/s)'].values
y0 = df0['AvgPowerUsage(W)'].values
y0_errors = df0['SDPowerUsage'].values
x0_errors = df0['SDAccessFrequency'].values

x1 = df1['AccessFrequency(line/s)'].values
y1 = df1['AvgPowerUsage(W)'].values
y1_errors = df1['SDPowerUsage'].values
x1_errors = df1['SDAccessFrequency'].values

x2 = df2['AccessFrequency(line/s)'].values
y2 = df2['AvgPowerUsage(W)'].values
y2_errors = df2['SDPowerUsage'].values
x2_errors = df2['SDAccessFrequency'].values

x3 = df3['AccessFrequency(line/s)'].values
y3 = df3['AvgPowerUsage(W)'].values
y3_errors = df3['SDPowerUsage'].values
x3_errors = df3['SDAccessFrequency'].values

# Scatter plot with error bars
plt.errorbar(x0, y0, yerr=y0_errors, fmt='o', color='blue', ecolor='blue', capsize=1, label="socket 0")
plt.errorbar(x1, y1, yerr=y1_errors, fmt='o', color='red', ecolor='red', capsize=1, label="socket 1")
plt.errorbar(x2, y2, yerr=y2_errors, fmt='o', color='green', ecolor='green', capsize=1, label="socket 2")
plt.errorbar(x3, y3, yerr=y3_errors, fmt='o', color='yellow', ecolor='yellow', capsize=1, label="socket 3")

# Least square
m0,c0 = np.polyfit(x0, y0, 1)
m1,c1 = np.polyfit(x1, y1, 1)
m2,c2 = np.polyfit(x2, y2, 1)
m3,c3 = np.polyfit(x3, y3, 1)
plt.plot(x0,m0*x0+c0, color='blue')
plt.plot(x1,m1*x1+c1, color='red')
plt.plot(x2,m2*x2+c2, color='green')
plt.plot(x3,m3*x3+c3, color='yellow')

# Display the plot
plt.title("ranker_coop1_socket_stream")
plt.xlabel("Access Frequency (line/s)")
plt.ylabel("Memory Power (W)")
plt.legend()
plt.grid(True)
plt.savefig("stream coop1 socket")
plt.show()
