# ENGR 5696
# JEFF EASTICK
print ("Assignment 1 Question 4")

import matplotlib.pyplot as plt
import numpy as np
import time

N = 50
n = 10

plt.ion()
canvas = np.zeros((N,N))

# Draw in an initial pattern with 1's
canvas[21,24] = 1
canvas[21,25] = 1
canvas[21,26] = 1


# Display matrix

def iteration(matrix):
    matrix2 = np.zeros((N,N))
    # cycle through matrix and check whether flips or not. create new matrix "canvas2" 
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            try:
                left = matrix[i-1,j]
                up = matrix[i,j-1]
                right = matrix[i,j+1]
                down = matrix[i+1,j]
            except IndexError as e:
                print("caught at border" + str(i) + str(j))

            total = left + up + right + down
            alive = total%2 # if even, alive = 0, make dead. if odd, alive = 1, make alive 
            matrix2[i,j] = 1 - alive
    return matrix2


plt.matshow(canvas)
plt.show()

for i in range(n):
    plt.matshow(canvas,fignum = False)
    plt.pause(0.001)
    canvas = iteration(canvas)
    plt.draw()




