# Example case from Svanberg MMA paper     
# Svanberg, "The method of moving asymptotes—a new method for structural optimization" (1993)

from nullspace_optimizer import EuclideanOptimizable, nlspace_solve
from nullspace_optimizer import symbolic_optimizable, bound_constraints_optimizable
        
C1 = 0.0624 
C2 = 1
    
@bound_constraints_optimizable(l=0)
@symbolic_optimizable
class P1(EuclideanOptimizable): 
    def x0(self):   
        return [5]*6    

    def J(self, x): 
        return C1*(x[1]+x[2]+x[3]+x[4]+x[5])    
        
    def H(self, x): 
        return [ 61/x[1]**3+37/x[2]**3+19/x[3]**3+7/x[4]**3+1/x[5]**3 - C2]
    
def solve_P1(): 
    problem = P1()  
    params = {'dt':1}
    results = nlspace_solve(problem,params)
    return results


if __name__=="__main__": 
    results = solve_P1()
    print("Optimal solution:")  
    print(f"x={results['x'][-1][1:]}")
