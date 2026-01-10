function [celasINV3D celas] = Compliance_ElasticityGIVEN(celas3D,typePROBLEM,StrainStressWith4Components)

 
%----------------------
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  G = E/2/(1+nu) ;
celasINV3D = inv(celas3D) ;
switch typePROBLEM
    case 'pstrain'
     %   celas3D = inv(celasINV3D)  ;
        if StrainStressWith4Components == 1
            rowcol = [1 2  6  3] ;  % Include also the normal stress along the z direction
        else
            rowcol = [1 2 6] ;
        end
        celas = celas3D(rowcol,rowcol) ;
    case 'pstress'
        rowcol = [1 2 6] ;
        celasINV = celasINV3D(rowcol,rowcol) ;
        celas = inv(celasINV) ;
    case '3D'
        celas = celas3D  ;
end