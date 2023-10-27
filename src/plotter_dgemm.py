import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Read CSV using pandas
df0 = pd.read_csv("../result/ranker_socket_dgemm-1,0_all.csv")
df1 = pd.read_csv("../result/ranker_socket_dgemm-0,7_all.csv")
df2 = pd.read_csv("../result/ranker_socket_dgemm-0,5_all.csv")
df3 = pd.read_csv("../result/ranker_socket_dgemm-0,4_all.csv")
df4 = pd.read_csv("../result/ranker_socket_dgemm-0,3_all.csv")
df5 = pd.read_csv("../result/ranker_socket_dgemm-0,2_all.csv")
    
# Extract values from dataframe
x0 = df0[' SocketId'].values
y0 = df0[' AvguJ/FLOP'].values
y0_errors = df0[' SDuJ/FLOP'].values

x1 = df1[' SocketId'].values
y1 = df1[' AvguJ/FLOP'].values
y1_errors = df1[' SDuJ/FLOP'].values

x2 = df2[' SocketId'].values
y2 = df2[' AvguJ/FLOP'].values
y2_errors = df2[' SDuJ/FLOP'].values

x3 = df3[' SocketId'].values
y3 = df3[' AvguJ/FLOP'].values
y3_errors = df3[' SDuJ/FLOP'].values

x4 = df4[' SocketId'].values
y4 = df4[' AvguJ/FLOP'].values
y4_errors = df4[' SDuJ/FLOP'].values

x5 = df5[' SocketId'].values
y5 = df5[' AvguJ/FLOP'].values
y5_errors = df5[' SDuJ/FLOP'].values

# Scatter plot with error bars
plt.errorbar(x0, y0, yerr=y0_errors, fmt='o-', color='blue', ecolor='blue', capsize=2, label="No power limit (225W)")
plt.errorbar(x1, y1, yerr=y1_errors, fmt='o-', color='red', ecolor='red', capsize=2, label="157 W")
plt.errorbar(x2, y2, yerr=y2_errors, fmt='o-', color='green', ecolor='green', capsize=2, label="112 W")
plt.errorbar(x3, y3, yerr=y3_errors, fmt='o-', color='yellow', ecolor='yellow', capsize=2, label="90 W")
plt.errorbar(x4, y4, yerr=y4_errors, fmt='o-', color='orange', ecolor='orange', capsize=2, label="67 W")
plt.errorbar(x5, y5, yerr=y5_errors, fmt='o-', color='pink', ecolor='pink', capsize=2, label="45 W")

plt.xticks(x0)

# Display the plot
plt.title("ranker_coop1_socket_dgemm")
plt.xlabel("Socket Id")
plt.ylabel("Avg uJ/FLOP")
plt.legend()
plt.grid(True)
plt.savefig("Dgemm coop1 socket")
plt.show()
