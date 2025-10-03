from nullspace_optimizer import nlspace_solve, EuclideanOptimizable 
from nullspace_optimizer import symbolic_optimizable

@symbolic_optimizable
class basicProblem(EuclideanOptimizable):
    def x0(self):
        return [1.5, 2.25]

    def J(self, x):
        return x[1]+0.3*x[0]

    def H(self, x):
        return [-x[1]+1.0/x[0], -(3-x[0]-x[1])]

def solve_basic_problem(**other_params):
    params = {'dt': 0.3}
    results = nlspace_solve(basicProblem(), params)
    return results

def main(**options):
    import nullspace_optimizer.examples.basic_examples.utils.draw as draw
    from nullspace_optimizer import undecorate

    results = solve_basic_problem(**options)
    print(f"Optimum : {results['x'][-1]}")

    draw.ion()
    draw.drawProblem(undecorate(basicProblem)(), XLIM=[0.2, 2.8], YLIM=[
                     0.2, 2.8], resolution=200)
    draw.drawData(results, 'Symbolic NLSPACE', 'blue')
    input("Press any key")
        
if __name__=="__main__":    
    main()
