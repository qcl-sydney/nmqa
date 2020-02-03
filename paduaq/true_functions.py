import numpy as np

def true_function(X, Y, d=None):
    '''Supports broadcasting; however riskanalysis.EngineeredTruth uses true_function
    for point values of X & Y so global operations for vector X and Y inside true function will fail. '''

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
        # [min, max] = 0 <= [0.8461838788653079, 3.1315926535897933] < pi over [-1, 1]^2
        result = np.pi*np.exp(-1.*((X+1)**2./5 + Y**2./2.)) -0.01
        return result

    if d == 'franke': # failed. Yields constant field under EngineeredTruth.

        # domain of franke function [0,1] --> [-1, 1]^2
        X = (X + 1.0)*0.5 
        Y = (Y + 1.0)*0.5
        
        # [min, max] = 0 <= [0.0008751540522927725 3.1395926535897933] < pi over [-1, 1]^2
        result =  0.75*np.exp(-((9*X - 2)**2 )/4 -((9*Y - 2)**2 )/4)
        result += 0.75*np.exp(-((9*X + 1)**2 )/49. -((9*Y + 1))/10)
        result +=  0.5*np.exp(-((9*X - 7)**2 )/4 -((9*Y - 3)**2 )/4)
        result +=  -0.2*np.exp(-((9*X - 4)**2 ) -((9*Y - 7)**2 ))
        
        result = np.pi * (result / np.max(result)) - 0.002 
        # returns constant value of 3.1 under EngineeredTruth.
        return result
        
    if d == 'franke_2':
        # domain of franke function [0,1] --> [-1, 1]^2
        X = (X + 1.0)*0.5 
        Y = (Y + 1.0)*0.5
        
        # [min, max] = 0 <= [0.0008751540522927725 3.1395926535897933] < pi over [-1, 1]^2
        result =  0.75*np.exp(-((9*X - 2)**2 )/4 -((9*Y - 2)**2 )/4)
        result += 0.75*np.exp(-((9*X + 1)**2 )/49. -((9*Y + 1))/10)
        result +=  0.5*np.exp(-((9*X - 7)**2 )/4 -((9*Y - 3)**2 )/4)
        result +=  -0.2*np.exp(-((9*X - 4)**2 ) -((9*Y - 7)**2 ))
        
        result = result * 2.5 + 0.0002 # local rescaling under EngineeredTruth.
        return result    

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
