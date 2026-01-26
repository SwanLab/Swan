import sys
sys.path.append("/home/joseantonio/Documentos/GitHub/Swan/Florian/NonLinearNull/03_NullSpace_HJ_Compl_VolPer_Iterative/")
import matplotlib.pyplot as plt

from fun03_NullSpace_HJ_CVP_Iterative import FunctionCase03

maxItj  = [1,2,4,8]
stepHJ  = 0.8
No      = 10
maxIter = 100

for i in range(len(maxItj)):
    FunctionCase03(maxItj[i],stepHJ,No,maxIter)
    plt.close('all')