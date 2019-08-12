import numpy as np

def true_function(X, Y, d=None):

    if d == 'cheb2fun':
        f_data = np.cos(np.exp(2*X+Y))*np.sin(Y)
        return (f_data + 1.)*np.pi/2.

    if d == 'lin':
        f_data = Y*0.5 + X*0.5 
        return (f_data + 1.)*np.pi/2.

    if d == 'exp':
        print "true_function not normalised to [0, pi]"
        return np.exp(-X) 

    if d == 'gss': 
        print "true_function not normalised to [0, pi]"
        return 2.0*np.exp((X-0.25)**2/5. + Y**2/30.)

    if d == 'sinc':
        print "true_function not normalised to [0, pi]"
        return np.sin(X**2 + Y**2) / (X**2 + Y**2)

    return np.pi * (X + Y + 0.1)/ (X + Y + 0.1)

def generate_data_qubits_coords(number, flag='random'):
    if flag == 'uniform':
        approxgrid = int(np.sqrt(number))
        data_x, data_y = np.meshgrid(np.linspace(-1, 1, approxgrid), np.linspace(-1, 1, approxgrid))
        data_qubits = zip(data_x.flatten(), data_y.flatten())
    if flag == 'random':
        data_qubits = np.random.rand(number, 2)*2.0 - 1.0 # inside unit square
    return data_qubits
