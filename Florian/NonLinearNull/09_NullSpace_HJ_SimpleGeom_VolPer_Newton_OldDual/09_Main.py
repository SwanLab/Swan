import sys
sys.path.append("/home/joseantonio/Documentos/GitHub/Swan/Florian/NonLinearNull/09_NullSpace_HJ_SimpleGeom_VolPer_Newton_OldDual/")
import matplotlib.pyplot as plt

from fun09_NullSpace_HJ_SimpleGeom_VP_Newton_OldDual import FunctionCase09

case    = ["Newton"]
No      = [25,25]
maxIter = 150

for i in range(len(case)):
    FunctionCase09(case[i],No[i],maxIter)
    plt.close('all')