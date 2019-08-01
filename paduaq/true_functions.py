import numpy as np

def true_function(X, Y, d=None):

    if d == 'cheb2fun':
        return np.cos(np.exp(2*X+Y))*np.sin(Y)

    if d == 'lin':
        return Y*0.5 + X*0.5

    if d == 'exp':
        return np.exp(X + 0.0*Y)

    if d == 'gss':
        return 2.0*np.exp((X-0.25)**2/5. + Y**2/30.)

    if d == 'sinc':
        return np.sin(X**2 + Y**2) / (X**2 + Y**2)

    return np.pi * (X + Y + 0.1)/ (X + Y + 0.1)

def generate_data_qubits_coords(number, flag='random'):
    if flag == 'uniform':
        approxgrid = int(np.sqrt(number))
        data_x, data_y = np.meshgrid(np.linspace(-1, 1, approxgrid), np.linspace(-1, 1, approxgrid))
        data_qubits = np.asarray(zip(data_x.flatten(), data_y.flatten()))
    if flag == 'random':
        data_qubits = np.random.rand(number, 2)*2.0 - 1.0 # inside unit square
    return data_qubits
