import numpy as np;

def mean(x):
    return np.sum(x)/np.size(x);

def std(x):
    m = mean(x);
    output = np.linalg.norm(m-x)/np.size(x);
    return output;



a = np.zeros(4);

for index, number in enumerate([1/100,1/80,1/60,1/40]):
    a[index] = 2.3-number;

print(a)
print(std(a))
print(mean(a))