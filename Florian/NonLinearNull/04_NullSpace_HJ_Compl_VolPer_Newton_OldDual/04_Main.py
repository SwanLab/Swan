import sys
sys.path.append("/home/joseantonio/Documentos/GitHub/Swan/Florian/NonLinearNull/04_NullSpace_HJ_Compl_VolPer_Newton_OldDual/")
import matplotlib.pyplot as plt

from fun04_NullSpace_HJ_CVP_Newton_OldDual import FunctionCase04

case    = "Original"
No      = 10
maxIter = 10

FunctionCase04(case,No,maxIter)
plt.close('all')