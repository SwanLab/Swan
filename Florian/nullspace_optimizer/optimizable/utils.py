import numpy as np  

def pack_alphas(alphas,p,q):
    if alphas is None:
        alphas = [1]*(p+q)
    alphas = np.asarray(alphas)
    if len(alphas) != (p+q):
        alphas = np.hstack((alphas,[1]*(p+q-len(alphas))))
    return alphas

def undecorate(instance):   
    if hasattr(instance,'_undecorated'):    
        return undecorate(instance._undecorated) 
    else:   
        return instance
