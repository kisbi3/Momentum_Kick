import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import multiprocessing
import time

# def function(x, y):
#     z = x*x + y
#     return z




# if __name__=="__main__":
#     x = np.arange(1, 6, 1)
#     y = np.arange(1, 6, 1)
#     print(x)
#     print(function(x, y))
#     pool = multiprocessing.Pool(processes=10)
#     res = pool.starmap(function, zip(range(1, 6), range(1, 6)))
#     print(res)


# x = cp.array([[0,1,2,3], [4,5,6,7]])
x = cp.array([0,1,2,3])
y = cp.array([[0,1,2,3], [5,6,7,8]])

print(cp.shape(x), cp.shape(y))
print(x, y)
print(x*y)