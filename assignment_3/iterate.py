import numpy as np 

def iterate_z(x,y,max_iter):
    """Iterates a complex number c through the equation z_{i} = z_{i+1}^2+ c a maximum
    number of times until the point diverges or not.
    
    Parameters:  
    x - float 
        the real part of the complex number z
    y - float 
        the imaginary part of the complex number z
    max_iter - int
        maximum number of iterations
    
    Returns: 
    c -  complex number 
         the original complex number 
    num - int
          the number of iterations need for the point to diverge. If 0, point did not diverge
    
    """
    #generate complex number 
    c = x + y*1j
    #our initial z
    z0 = 0
    #iterating the point 
    for i in range(max_iter):
        z0 = z0**2 + c 
        #Divergence condition 
        if abs(z0) > 10**(100):
            #returns the complex point and iteration 
            return c,i
    return c, 0
