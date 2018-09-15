# ENGR 5696
# JEFF EASTICK
print ("Assignment 1 Question 4")

import matplotlib.pyplot as plt
import numpy as np

N = 50
n = 40

plt.ion()
canvas = np.zeros((N,N))

# Draw an initial pattern with 1's
canvas[20,20] = 1
canvas[20,21] = 1
canvas[20,22] = 1
canvas[20,23] = 1
canvas[20,24] = 1



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



# Main code
plt.matshow(canvas)
plt.show()

for i in range(n):
    plt.matshow(canvas,fignum = False)
    plt.pause(0.01)
    canvas = iteration(canvas)




