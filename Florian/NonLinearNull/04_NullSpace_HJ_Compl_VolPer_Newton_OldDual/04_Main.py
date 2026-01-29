import sys
sys.path.append("/home/joseantonio/Documentos/GitHub/Swan/Florian/NonLinearNull/04_NullSpace_HJ_Compl_VolPer_Newton_OldDual/")
import matplotlib.pyplot as plt

from fun04_NullSpace_HJ_CVP_Newton_OldDual import FunctionCase04

case    = ["Original","Newton"]
No      = [40,40]
maxIter = 100

for i in range(len(case)):
    FunctionCase04(case[i],No[i],maxIter)
    plt.close('all')