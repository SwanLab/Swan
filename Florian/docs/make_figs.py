from utils import savefig, exec1
from nullspace_optimizer.examples.basic_examples import ex02_basic_problem_3, ex05_multiple_shootings, ex06_multiple_shootings2
import matplotlib.pyplot as plt
import nullspace_optimizer.examples.basic_examples.utils.draw as draw

results = ex02_basic_problem_3.solve_problem_parab()
draw.drawProblem(ex02_basic_problem_3.parabProblem(), XLIM=[-3.5, 5], YLIM=[-1.3, 3.2])
draw.drawData(results, 'NLSPACE', 'blue', loc='lower right')
draw.remove_legend()
savefig('source/img/ex02.png', density=200)
plt.close('all')

results = ex05_multiple_shootings.run_problems()



draw.drawProblem(ex05_multiple_shootings.problemeSimple([0, 0]), [-1.5, 1.5], [-1.5, 1.5])
for i, r in enumerate(results):
    draw.drawData(r, f'x{i+1}', f'C{i}', x0=True, xfinal=True, initlabel=None)
draw.remove_legend()
        
savefig('source/img/ex05.png', density=200)
plt.close('all')
    
results = ex06_multiple_shootings2.run_problems()

draw.drawProblem(ex06_multiple_shootings2.problemeSimple2([0, 0]), [-1.5, 1.5], [-1.5, 1.5])
for i, r in enumerate(results):
    draw.drawData(r, None, f'C{i}', x0=True, xfinal=True)
draw.remove_legend()
savefig('source/img/ex06.png', density=200)
    
plt.close('all')
    

from nullspace_optimizer import EuclideanOptimizable, nlspace_solve
class problemSimple(EuclideanOptimizable):
   # Initialization
   def x0(self):
       return [1.25,0]

   # Objective function
   def J(self, x):
       return (x[0]+1)**2+(x[1]+1)**2

   # Inequality constraint
   def H(self, x):
       return [x[0]**2+x[1]**2-1**2,
               x[0]+x[1]-1,
               -x[1]-0.7]

   # Row vector with the gradient of J
   def dJ(self, x):
       return [2*(x[0]+1), 2*(x[1]+1)]

   # Jacobian matrix of H
   def dH(self, x):
       return [[2*x[0], 2*x[1]],
               [1, 1],
               [0, -1]]

params = {'dt': 0.1}
opt_results = nlspace_solve(problemSimple(), params)

draw.drawProblem(problemSimple(), [-1.5, 1.5], [-1.5, 1.5])
draw.drawData(opt_results, f'x{0}', f'C{0}', x0=True, xfinal=True, initlabel=None)
        
savefig('source/img/ex_doc.png', density=200)
plt.close('all')

exec1('mogrify -trim source/img/*.png')
