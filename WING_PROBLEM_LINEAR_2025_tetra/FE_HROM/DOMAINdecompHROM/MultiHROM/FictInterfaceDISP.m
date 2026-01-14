function  [Vall,MESH,DATA,INFOPLOTMESHBND,FluctuationFacesGLO] = FictInterfaceDISP(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline)
% Displacement  of fictitious interfaces
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/README_HOMOG.mlx
% JAHO 12-FEB-2023
%5------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

 
INFOPLOTMESHBND = [] ; 
FluctuationFacesGLO = [] ; 
% ---------------------------------------------------------------------------------------------
DATAcommon = DefaultField(DATAcommon,'TypeFunctionDisplacementInterfaces','QUADRILATERAL_LINEAR') ;
 switch DATAcommon.TypeFunctionDisplacementInterfaces
     case ''
    case 'QUADRILATERAL_LINEAR'
    [Vall,Mintf,MESH,DATA,FluctuationFacesGLO] = FictInterfaceDISP_QUADlin(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline) ; 
    
     case 'QUADRILATERAL_QUADRATIC'
    [Vall,Mintf,MESH,DATA,INFOPLOTMESHBND] = FictInterfaceDISP_QUADq(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline) ; 
    
    case 'PLATE_QUAD_LINEAR'
     [Vall,Mintf,MESH,DATA] = FictInterfaceDISP_PLATEquadLIN(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline) ; 
      case 'TRIANGULAR_LINEAR'
        [Vall,Mintf,MESH,DATA] = FictInterfaceDISP_TRIlin(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline) ; 
       case {'HEXAHEDRA_LINEAR','HEXAHEDRA_8nodes'}  % Alternative names
    [Vall,Mintf,MESH,DATA] = FictInterfaceDISP_HEXAlin(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline) ;    
        
 case {'HEXAHEDRA_QUADRATIC','HEXAHEDRA_27nodes'}  % Alternative names
   %  disp('Rather use ')
    [Vall,MESH,DATA,INFOPLOTMESHBND] = FictInterfaceDISP_HEXAquad(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline) ;    
    case 'HEXAHEDRA_20nodes'
    [Vall,Mintf,MESH,DATA] = FictInterfaceDISP_HEXA20NODES(PhiDEF,PsiDEFf,MESH,DATAcommon,DATA,Mintf,PhiRB,PsiRBf,DATAoffline) ;    

    
    
    otherwise
        error('Option not implemented yet')
end

