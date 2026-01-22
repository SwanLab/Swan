import sys
sys.path.append("/home/joseantonio/Documentos/GitHub/Swan/Florian/NonLinearNull/03_NullSpace_HJ_Compl_VolPer_Iterative/")
import matplotlib.pyplot as plt

from fun03_NullSpace_HJ_CVP_Iterative import FunctionCase03

maxItj  = 1
stepHJ  = 0.3
No      = 20
maxIter = 100

FunctionCase03(maxItj,stepHJ,No,maxIter)
plt.close('all')