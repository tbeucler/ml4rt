""" 
S. Shamekh, ss6287@columbia.edu
11/10/2021


This script contains functions for computing
Rieman Liouville(RL) fractional integral and   
Grunwald (GL) fractional integrroderivative.

RL functions: take a an interval as well as alpha and 
return the fractional integral of all points. To compute the 
derivaties one can deirivate the integrals. 
More details:  Numerical Approaches to Fractional Integrals and Derivatives: A Review

GL functions: take un interval and return the fractional integoderivative for one side of the interval, 
it has to be combined with an outside loop to get the fractional integroderivatives for all points.
More details:  The fractional calculus, K. B. Oldham, J. Spanier; 1974

"""


def RL_left_integral(data,dx,alpha):
#     """ Computes the Reiman-Liouville fractional integral of a function at a point.
#        based on equation 20 and 21 of https://www.mdpi.com/2227-7390/8/1/43
#        Parameters
#        ==========
#         data : 2D array (time,height) to be differintegrated.
       
#         alpha : float and positive
#             The order of the differintegral to be computed. Height is (a,b)
#         dx   :  the descretization
        
#         return : left-sided fractional integral for all points (a,b)
        

#     """
    fractional = np.zeros_like(data)
    
    for j in range(1,data.shape[1]):
        
        coefficient = np.zeros((j+1))
        N = j+1
        coefficient[0] = np.power(j-1,alpha+1) - (j-1-alpha)*np.power(j,alpha)
        coefficient[N-1] = 1

        for k in range(1,N-1):
            coefficient[k] = np.power(j-k+1,alpha+1) - 2*np.power(j-k,alpha+1) + np.power(j-k-1,alpha+1)
            
        fractional[:,j] = np.power(dx,alpha)*np.dot(data[:,:j+1],coefficient)/gamma(alpha+2)
    return fractional


def rRL_right_integral(data,dx,alpha):
#      """ Computes the Reiman-Liuville fractional integral of a function at a point.
       
#        Parameters
#        ==========
#         data : 2D array (time,height) to be differintegrated. height is (a,b)
       
#         alpha : float and positive 
#             The order of the differintegral to be computed.
#         dx   :  the descretization
        
#         return : right-sided fractional integral at all points (a,b)

#     """

    
    fractional = np.zeros_like(data)
    for j in range(0,data.shape[1]-1):
              
        L = len(data[0,j:])
        N = L-1
        coefficient = np.zeros((L))
        coefficient[0] = 1
        coefficient[L-1] =  np.power(N-1,alpha+1) - (N-1-alpha)*np.power(N,alpha)

        for k in range(1,L-1):
            coefficient[k] = np.power(k+1,alpha+1) - 2*np.power(k,alpha+1) + np.power(k-1,alpha+1)
        fractional[:,j] = np.power(dx,alpha)*np.dot(data[:,j:],coefficient)/gamma(alpha+2)
    return fractional




def GL_right(data,alpha,dx):
    
#       """ Computes the Grunwald fractional integral/derivative of a function at a point.
#        based on Meerschaert and Tadjeran 2006, eq 1.6-1.8
       
#        Parameters
#        ==========
#         data : 2D array (time,height) to be differintegrated.  height is (a,z)

#         alpha : float; positive--> derivative , negative--> integral 
#             The order of the differintegral to be computed.
#         dx   :  the descretization
        
#         return : fractional integral/derivative at z

#     """
    length = data.shape[1]
    coeff2 = np.ones((length))
    
    for k in range(0,length):
        if k == 0 :
            coeff2[k] = 1
        else : 
            for i in range(1,k+1):
                coeff2[k] *= (alpha-(i-1))/(i)
            coeff2[k]*=np.power(-1,k)
        
    return np.dot(np.flip(data,axis=1),coeff2)/(dx**alpha)

def GL_left(data,alpha,dx):
    
#       """ Computes the Grunwald fractional integral/derivative of a function at a point.
#        based on Meerschaert and Tadjeran 2006, eq 1.6-1.8
       
#        Parameters
#        ==========

#         data : 2D array (time,height) to be differintegrated. height is (z,b)
        
#         alpha : float; positive--> derivative , negative--> integral 
#         The order of the differintegral to be computed.
        
#         dx   :  the descretization
        
#         return : fractional integral/derivative at z

#     """
    length =data.shape[1]
    coeff2 = np.ones((length))
    
    for k in range(0,length):
        if k == 0 :
            coeff2[k] = 1
        else:
            for i in range(1,k+1):
                coeff2[k] *= (alpha-(i-1))/(i)
            coeff2[k]*=np.power(-1,k)
    return np.dot(data,coeff2)/(dx**alpha)

