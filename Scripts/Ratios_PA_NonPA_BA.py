#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 15:38:24 2024

@author: laura
"""
import numpy as np
import matplotlib.pyplot as plt
#%%

path_fig = '/Users/laura/Masterthesis/Figures/'
#%%

counts_pa = [0, 13, 1, 0, 0, 0, 77, 292, 639, 21, 98, 148, 1, 0, 0, 7, 0, 36, 280, 530, 155, 0, 0, 0, 38, 0, 2, 38, 0, 146, 31, 3633, 1597, 2, 0, 0, 0, 20, 29, 1, 8, 1, 590, 441, 25, 94, 281, 44, 81, 3, 15, 11, 0, 22, 1499, 3182, 293, 235, 161, 226, 7, 17, 0, 0, 15, 61, 74, 706, 263, 0, 0, 136, 17, 0, 26, 0, 124, 0, 0, 76, 45, 35, 631, 9, 1, 17, 32, 39, 0, 0, 0, 24, 12, 77, 89, 9, 10, 73, 501, 64, 3, 0, 0, 11, 199, 44, 0, 0, 0, 0, 27, 9, 0, 0, 8, 1883, 193, 2, 76, 0, 61, 52, 19, 24, 11, 0, 25, 74, 91, 451, 0, 0, 0, 139, 273, 28, 0, 0, 764, 7, 166, 2, 0, 0, 0, 0, 0, 5, 0, 0, 276, 561, 527, 203, 31, 142, 0, 0, 9, 0, 0, 136, 0, 274, 16, 0, 0, 14, 303, 10, 130, 129, 58, 8, 238, 487, 65, 31, 0, 0, 0, 10, 0, 0, 0, 3, 12, 4330, 237, 48, 9, 128, 160, 15, 83, 145, 9, 200, 114, 831, 195, 2872, 10, 0, 31, 5, 0, 0, 0, 0, 0, 984, 72, 58, 0, 0, 835, 3, 107, 8, 0, 0, 0, 2, 26, 0, 0, 6, 66, 63, 0, 0, 0, 0, 319, 34, 155, 0, 77, 79]
counts_nonpa = [0, 0, 0, 0, 0, 0, 52, 366, 2302, 157, 18, 845, 16, 0, 0, 10, 0, 69, 1197, 1886, 304, 0, 0, 0, 1, 0, 10, 1, 0, 311, 3175, 21628, 852, 46, 0, 0, 0, 1, 1, 176, 2, 45, 1117, 612, 162, 86, 13, 0, 0, 72, 36, 16, 0, 488, 4301, 13013, 1941, 653, 4, 64, 5, 0, 0, 0, 14, 343, 143, 1946, 160, 0, 0, 0, 0, 0, 8, 0, 2, 0, 31, 235, 21, 30, 835, 54, 0, 1, 33, 11, 0, 0, 69, 56, 10, 89, 2, 86, 9, 31, 429, 29, 8, 12, 123, 97, 902, 344, 0, 0, 0, 0, 0, 82, 0, 0, 342, 3016, 677, 63, 0, 0, 78, 20, 7, 1, 0, 85, 329, 250, 94, 1260, 2, 0, 1, 50, 342, 246, 0, 0, 592, 80, 1002, 10, 0, 0, 0, 0, 0, 0, 14, 48, 428, 1474, 1245, 280, 8, 109, 0, 0, 13, 0, 3, 166, 27, 156, 125, 0, 6, 151, 5, 0, 63, 341, 289, 70, 373, 1288, 103, 67, 0, 0, 0, 7, 0, 7, 0, 0, 134, 6161, 1151, 140, 62, 2, 52, 10, 105, 67, 40, 6774, 3758, 5977, 2363, 19694, 7, 0, 6, 0, 0, 0, 0, 0, 0, 399, 29, 91, 0, 90, 652, 0, 93, 27, 7, 0, 1544, 70, 415, 8, 0, 2, 14, 7, 0, 0, 0, 144, 1213, 717, 3118, 3, 84, 27]

pixarea = 0.0625
0.00000625 * 280

list1 = sum(counts_pa)
list2 = sum(counts_nonpa)
allpa = list1 + list2

list1/allpa
list2/allpa


ratios = []
for i in range(len(counts_pa)):
    if counts_nonpa[i] != 0:
        ratios.append(counts_pa[i] / counts_nonpa[i])
    else:
        ratios.append(float('inf'))  # or handle division by zero as you prefer

above_one_count = 0
below_one_count = 0
inf_count = 0

for ratio in ratios:
    if ratio == float('inf'):
        inf_count += 1
    elif ratio > 1:
        above_one_count += 1
    elif ratio < 1:
        below_one_count += 1
        
filtered_ratios = [ratio for ratio in ratios if ratio != float('inf')]
if filtered_ratios:
    max_ratio = max(filtered_ratios)
else:
    max_ratio = None 
    
#%%

yearly_sum_list1 = np.sum(np.reshape(counts_pa, (20, 12)), axis=1)
yearly_sum_list2 = np.sum(np.reshape(counts_nonpa, (20, 12)), axis=1)

# Step 2: Create a list of years from 2001 to 2020
years = list(range(2001, 2021))

#%%
# Step 3: Plot both time series on the same plot
plt.figure(figsize=(10, 6))
plt.plot(years, yearly_sum_list1, label='Forest in PA', marker='o', color='darkgreen')
plt.plot(years, yearly_sum_list2, label='Forest in Non-PA', marker='o', color='orange')

# Step 4: Customize the plot
plt.xlabel('Year')
plt.ylabel('Burned Pixels Count')
plt.title('Yearly Burned Pixels (2001-2020) in Portuguese forests')
plt.legend()
plt.grid(True)
plt.xticks(years, rotation=45)  # Ensure all years are visible on the x-axis
plt.tight_layout()

plt.savefig(path_fig+'BA_count_forest_years.png', dpi=250)









