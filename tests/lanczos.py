#Testing script for Lanczos Algorithm
import numpy as np

N = 9
#Create simple symmetric matrix
A = np.zeros((N, N))
for i in range(N):
    for j in range(N):
        A[i, j] = 1/(i + j + 1)
    A[i, i] = 1/(i + 1)

X0 = np.ones(N)
X0 = X0/np.linalg.norm(X0)
def lanczos(A, X0):
    V = np.zeros((N, N))

    def matvec(v):
        return A.dot(v)
    
    def mgs(index):
        result = V[:, index]
        for i in range(index):
            result = result - np.dot(result, V[:, i]) * V[:, i]
        result = result/np.linalg.norm(result)
        return result

    V[:, 0] = X0
    V[:, 0] = V[:, 0]/np.linalg.norm(V[:, 0])
    alphas  = np.zeros(N)
    betas   = np.zeros(N)
    for i in range(N):
        if (i > 0 and i % 2 == 0):
            V[:, i-1] = mgs(i-1)
            V[:, i] = mgs(i)
        
        v = matvec(V[:, i])
        alphas[i] = np.dot(V[:, i], v)        
        if (i == 0):
            v = v - alphas[i] * V[:, i]
        else:
            v = v - alphas[i] * V[:, i] - betas[i-1] * V[:, i-1]
        betas[i] = np.linalg.norm(v)
        if (i < N - 1):
            V[:, i+1] = v/betas[i]

    #Construct resultant tridiagonal matrix
    T = np.zeros((N, N))
    for i in range(N):
        T[i, i] = alphas[i]
        if (i < N - 1):
            T[i, i+1] = betas[i]
            T[i+1, i] = betas[i]
    return T, alphas, betas

T, alphas, betas = lanczos(A, X0)
print(np.linalg.eigh(A)[0])
print(np.linalg.eigh(T)[0])
#print(alphas)
#print(betas)
print(np.linalg.eigh(A)[1][:, -3:])
print((A @ np.linalg.eigh(A)[1][:, -1])/np.linalg.eigh(A)[1][:, -1])
#X = np.array([0.679183, 0.456510, 0.332463, 0.263441, 0.219243, 0.188302, 0.165318, 0.147511, 0.133276])
#print((A @ X) / X)



#Random A matrix
""" A = np.random.rand(N, N)
#force symmetry
for i in range(N):
    for j in range(i):
        A[i, j] = A[j, i]

X = np.random.rand(N,3)

print(X.T @ A @ X)
 """