function [VAR,celastST,FgradST,detFgrad] = StressesFromDisplacementsCABLES(OPERFE,VAR,MATPRO,DATA,VARint_n)


% See  /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/FIBREGY_PROJECT_2022/04_MOORING_PROBLEMS/01_STATICgrav.mlx
if nargin == 0
    load('tmp.mat')
end
detFgrad = [] ;
% Equivalent deformation gradient (F)
ndim = DATA.MESH.ndim ;
ndof = size(OPERFE.Bst,2) ;
FgradST =  OPERFE.Bst*VAR.DISP(1:ndof,:)  ; %+ OPERFE.IDENTITY_F ;
idim = 1 ;
LOCROWS = idim:ndim:size(FgradST,1) ;
FgradST(LOCROWS,:) = 1+FgradST(LOCROWS,:) ;

% 3. Green-Lagrante strains at all Gauss points (equivalent for cables)
%VAR.GLSTRAINS = StrainGreenLagrange(FgradST,DATA.MESH.ndim) ;

FtF = zeros(size(VAR.STRAIN)) ;
for idim = 1:ndim
    FtF = FtF + FgradST(idim:ndim:end,:).^2 ;
end
 
VAR.STRAIN =sqrt(FtF)-1; ;


% Tension = Area*Stress

[stress,der_stress] = StressesFromStrainCABLE(MATPRO,VAR.STRAIN,DATA) ; 
VAR.TENSION = MATPRO.AREA.*stress;

T = VAR.TENSION./(1+VAR.STRAIN) ; %

for idim = 1:ndim
    VAR.TENSIONV(idim:ndim:end,:) = T.*FgradST(idim:ndim:end,:); 
end
celastST = []; 
if  DATA.CALC_CTANG == 1
    
    % Constitutive tangent matrix
    % Geometric contribution
    celastST = zeros(size(VAR.TENSIONV,1),ndim) ;
    
    
    if DATA.AVOID_NEGATIVE_GEOMETRIC_MATRICES ==1
        IndCompression =  find(T<=0)  ;
        T(T<=0) =  MATPRO.InitialTension(IndCompression) ;
    end
    
    for idim = 1:ndim
        celastST(idim:ndim:end,idim) = T;
    end
    
    % Material contribution
    % ---------------------
    Ebar = der_stress + (der_stress.*VAR.STRAIN - stress) ; 
    factor_E = (MATPRO.AREA.*Ebar)./((1+VAR.STRAIN).^3) ;
    for idim = 1:ndim
        for jdim = 1:ndim
            celastST(idim:ndim:end,jdim) = celastST(idim:ndim:end,jdim) + factor_E.*(FgradST(idim:ndim:end).* FgradST(jdim:ndim:end)) ;
        end
    end
    
end

 