# This file is part of nullspace_optimizer. 
#
# It is a modification of the 2018-2019 version initially developed with the    
# copyright of CNRS, Ecole Polytechnique and Safran.
#
# This version is maintained by the research team of Prof. Florian Feppon 2023
# https://people.cs.kuleuven.be/~florian.feppon/software.html
#
# nullspace_optimizer is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# nullspace_optimizer is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is included below.
# For further information, see <http://www.gnu.org/licenses/>.
        
import numpy as np
import scipy.sparse as sp
from .optimizable import Optimizable
from .inout import display, colored
import hashlib

def hash(data): 
    __hash__ = hashlib.sha256() 
    if hasattr(data,'tobytes'):    
        __hash__.update(data.tobytes())
    else:
        __hash__.update(repr(data).encode())
    return __hash__.digest()

def finiteDiffCheck(problem,x,eps=1e-6,dx=None,seed=None):
    J = problem.J(x)
    dJ = problem.dJ(x)
    if seed:
        np.random.seed(seed)
    if not dx is None:
        dx = np.random.rand(len(dJ))
    dx = dx/(1e-15+np.linalg.norm(dx,np.inf))*eps
    newx = problem.retract(x, dx)
    newJ = problem.J(newx)
    checkJ = abs(newJ - J -dJ.dot(dx))/eps
    display(f"J={J}",color="magenta")
    display(f"newJ={newJ}",color="magenta")
    display("Numerical sensitivity:",color="magenta")
    display(f"(J(x+dx)-J(x)-dJ(dx))/||dx||={checkJ}",color="magenta")
        
def memoize(hash_function = None, store=5, func_name = None, debug=0):  
    """ 
    A decorator for memoization of functions. This allows to use    
    a saved value when calling several time a function on the same  
    argument.   
        
    **Usage**:  
        
    .. code:: python    
        
       @memoize() 
       def func(x): 
            return x+3  
            
       func(3) # Compute the function by calling func   
       func(3) # Use the already computed value 
       func(5) # Compute the function by calling func
        
    ``@memoize()`` can also be used with functions supporting multiple arguments.
    
    .. code:: python    
        
       @memoize() 
       def func(x1,x2): 
            return x1+x2  
        
    :param hash_function: a custom hash function for hashing function arguments.    
                          By default, the hash function is 
                          `hashlib.sha256 <https://docs.python.org/3/library/hashlib.html>`_    
                          applied to the either the string representation of the argument,  
                          or the bytes representation if the argument is a numpy array.
    :param store: Number of stored values in the hash table. 
    :type store: int
        
    :param func_name: the name of the function printed when debugging 
    :type func_name: str    
        
    :param debug: an integer tuning the level of verbosity. When greater than 0,    
                  information about whether the original function or its memoization    
                  is displayed.     
                    
    :type debug: int
    """
    if hash_function is None:   
        hash_function = hash
    save = dict()
    queue = []
    def memoize_impl(func):
        def newfunc(*args): 
            hashed = tuple(hash_function(arg) for arg in args)
            if not hashed in save:
                if debug>=1 and func_name:
                    message = colored("@memoize: ",color="magenta")+colored(func_name,attr="bold")+\
                        colored(" calling the function !",color="deep_pink_4c", attr="bold")   
                    if debug >=3:
                        message +=" (x="+",".join(map(str,args))+")"
                    print(message)
                save[hashed] = func(*args)
                queue.append(hashed)
                if len(queue)>store:    
                    save.pop(queue[0])  
                    queue.pop(0)
                if debug >= 3:
                    print(colored("@memoize: ",color="magenta")+f"{queue=}")

            elif debug>=2 and func_name:   
                message = colored("@memoize: ",color="magenta")+colored(func_name,attr="bold")  \
                    +colored(" using stored value.",color="magenta")    
                if debug >=3:
                    message += " (x="+",".join(map(str,args))+")"
                print(message)
            return save[hashed]
        return newfunc
    return memoize_impl

class OptimizationState():
    """ 
    An optimization state is working like a dictionary object, but  
    each element can be accessed only once.
    """

    def __init__(self, update_function, hash_function = None):
        self.__entries__ = dict()
        self.__update_function__ = update_function
        self.__lastx__ = None
        if hash_function is None:   
            self.hash = hash
        
    def __call__(self, x, key):
        if not self.hash(x) == self.__lastx__:
            self.__entries__ = self.__update_function__(x)
            self.__lastx__ = self.hash(x)
        return self.__entries__[key]
        
def eval_problem(problem : Optimizable, x):    
    """ 
    Shortcut function to obtain objective, constraints, and sensitivities.  
        
    :param problem: an :class:`Optimizable` object  
    :param x: an input point on which to evaluate the problem   
    :return: ``J, G, H, dJ, dG, dG``: objective function, constraints, and sensitivities    

    .. note::   
        
       These values are not processed into numpy arrays or proper sparse formats.   
       See :func:`optimizers.nullspace.utils.pack_constraints` for doing so.
    """
    J = problem.G(x)
    G = problem.G(x)
    H = problem.H(x)   
    dJ = problem.dJ(x)
    dG = problem.dG(x)     
    dH = problem.dH(x) 
    return J, G, H, dJ, dG, dH
