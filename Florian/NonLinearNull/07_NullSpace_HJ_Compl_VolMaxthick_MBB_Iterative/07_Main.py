import sys
sys.path.append("/home/joseantonio/Documentos/GitHub/Swan/Florian/NonLinearNull/07_NullSpace_HJ_Compl_VolMaxthick_MBB_Iterative/")
import matplotlib.pyplot as plt

from fun07_NullSpace_HJ_CVT_Iterative import FunctionCase07

maxItj  = [1,5]
stepHJ  = 1
No      = 10
maxIter = 40

for i in range(len(maxItj)):
    FunctionCase07(maxItj[i],stepHJ,No,maxIter)
    plt.close('all')