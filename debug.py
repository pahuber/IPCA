# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 18:07:14 2019

@author: philipp
"""

import matplotlib.pyplot as plt

x = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80]
y_inner = [6.61, 5.39, 7.41, 7.75, 6.53, 5.83, 5.51, 5.35, 5.31, 5.26, 5.19, 5.12, 5.05, 5.02, 5.03]
y_outer = [10.03, 11.65, 10.85, 10.45, 9.85, 9.59, 9.80, 9.89, 9.96, 10.04, 10.10, 10.20, 10.25, 10.38, 10.51]


plt.plot(x, y_inner, color="b", marker=".", linestyle=":", label="inner")
plt.plot(x, y_outer, color="r", marker=".", linestyle=":", label="outer")

plt.title("SNR vs. Rank")
plt.xlabel("Rank")
plt.ylabel("SNR")
plt.legend()
plt.grid(axis="y")
#plt.show()
plt.savefig("snr.png", dpi=1000)