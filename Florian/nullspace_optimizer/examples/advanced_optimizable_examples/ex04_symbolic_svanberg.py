# Example case from Svanberg MMA paper     
# Svanberg, "The method of moving asymptotes—a new method for structural optimization" (1993)

from nullspace_optimizer import EuclideanOptimizable, nlspace_solve
from nullspace_optimizer import symbolic_optimizable, bound_constraints_optimizable
from sympy import sqrt
        
C1 = 1 
C2 = 0.124
    
@bound_constraints_optimizable(l=(None,0.2,0.1),u=(None,4.0,1.6))
@symbolic_optimizable
class P1(EuclideanOptimizable): 
    def x0(self):   
        return [0,1.5,0.5]

    def J(self, x): 
        return C1*x[1]*sqrt(1+x[2]**2)  
        
    def H(self, x): 
        return [ C2*sqrt(1+x[2]**2)*(8/x[1]+1/(x[1]*x[2]))-1,   
                 C2*sqrt(1+x[2]**2)*(8/x[1]-1/(x[1]*x[2]))-1]

def solve_P1(): 
    problem = P1()  
    params = {'dt':0.1}
    results = nlspace_solve(problem,params)
    return results


if __name__=="__main__": 
    results = solve_P1()
    print("Optimal solution:")  
    print(f"x={results['x'][-1][1:]}")

