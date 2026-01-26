import sys
sys.path.append("/home/joseantonio/Documentos/GitHub/Swan/Florian/NonLinearNull/03_NullSpace_HJ_Compl_VolPer_Iterative/")
import matplotlib.pyplot as plt

from fun03_NullSpace_HJ_CVP_Iterative import FunctionCase03

maxItj  = [1,2,4,8]
stepHJ  = 1
No      = 12
maxIter = 40

for i in range(len(maxItj)):
    FunctionCase03(maxItj[i],stepHJ,No,maxIter)
    plt.close('all')