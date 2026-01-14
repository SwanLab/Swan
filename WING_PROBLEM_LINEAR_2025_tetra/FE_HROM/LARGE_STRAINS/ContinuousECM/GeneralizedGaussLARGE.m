function   [xNEW,wNEW,DATAIN,ELEMENTS_xNEW,HISTORY_ITERATIONS,VAR_SMOOTH_FE,Ninterpolation] =  GeneralizedGaussLARGE...
    (W,MESH,DATAIN,HYPERREDUCED_VARIABLES,Nst,DATALOC,VAR_SMOOTH_FE)  ;


if nargin == 0
    load('tmp1.mat')
end

% xy --> Coordinates GAUSS points
COORg = MESH.COOR' ;
xGAUSSall = Nst*COORg(:) ;
ndim = size(MESH.COOR,2) ;
COORg = reshape(xGAUSSall,ndim,[])' ;
VAR_SMOOTH_FE.COORg  = COORg ; 

% Indexes integration points
z = HYPERREDUCED_VARIABLES.setPoints;
% REduced weights
w = HYPERREDUCED_VARIABLES.WdomRED ;
coorECM  = COORg(z,:) ;   
 % Exact Integral 
DATALOC.ExactIntegral =  HYPERREDUCED_VARIABLES.PHI'*W ; 

VAR_SMOOTH_FE.BasisIntegrand = HYPERREDUCED_VARIABLES.PHI ; 
VAR_SMOOTH_FE.setPoints = z ; 


DATALOC =DefaultField(DATALOC,'CONTROLLED_POINTS_WEIGHTS_APPROACH',0) ; 

if DATALOC.CONTROLLED_POINTS_WEIGHTS_APPROACH == 0
    % Approach before 2nd-January-2022
[xNEW,wNEW,DATALOC,ELEMENTS_xNEW,HISTORY_ITERATIONS,VAR_SMOOTH_FE,POLYINFO] ...
    = GenGauss2D3Dfebased_large(COORg,w,DATALOC,coorECM,VAR_SMOOTH_FE) ;

else
    % Approach after 2nd-January-2022
    % See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/
    % ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx
    
    [xNEW,wNEW,DATALOC,ELEMENTS_xNEW,HISTORY_ITERATIONS,VAR_SMOOTH_FE,POLYINFO] ...
    = GenGauss2D3DfebasedCONTROL(COORg,w,DATALOC,coorECM,VAR_SMOOTH_FE) ;
    
end


VAR_SMOOTH_FE.IndexECMini = z ; 
%VAR_SMOOTH_FE.Bhred_interp =  BdomRED_interp ; 


VAR_SMOOTH_FE.ELEMENTS_xNEW = ELEMENTS_xNEW ;
VAR_SMOOTH_FE.COORiniECM  = HISTORY_ITERATIONS.POINTS{1}; % iNITIAL SET OF ECM POINTS  
VAR_SMOOTH_FE.WEIGHTSiniECM  = HISTORY_ITERATIONS.WEIGHTS{1}; % iNITIAL SET OF ECM POINTS  
VAR_SMOOTH_FE.ELEMENTSiniECM  = HISTORY_ITERATIONS.ELEMENTS_CONTAINING_POINTS{1}; % iNITIAL SET OF ECM POINTS 

% In order to perform interpolations of any Gauss variable, it is necessary
% to have at one's dispossal the coefficients of the shape functions as
% well as the scaling factos of ELEMENTS_xNEW

%INTERPOLATION_INFO.COEFFSpolynomial = POLYINFO.COEFFSpolynomial(ELEMENTS_xNEW);
%INTERPOLATION_INFO.SCALING_VARIABLES.LENGTH = POLYINFO.SCALING_VARIABLES.LENGTH(ELEMENTS_xNEW);
%INTERPOLATION_INFO.SCALING_VARIABLES.LENGTH = POLYINFO.SCALING_VARIABLES.LENGTH(ELEMENTS_xNEW);
%INTERPOLATION_INFO.SCALING_VARIABLES.coorREF = POLYINFO.SCALING_VARIABLES.coorREF(ELEMENTS_xNEW,:);

% INTERPOLATION FUNCTIONS 
Ninterpolation = cell(length(ELEMENTS_xNEW),1) ; % 
for ipoint = 1:length(ELEMENTS_xNEW)
    % Loop over elements 
    % Coordinate of the point under consideration 
    xLOC = xNEW(ipoint,:) ; 
    ielemLOC = ELEMENTS_xNEW(ipoint) ; 
    % Length scaling 
    LelemSCALING = POLYINFO.SCALING_VARIABLES.LENGTH(ielemLOC,:) ; 
    % Reference point element 
    coorREF = POLYINFO.SCALING_VARIABLES.coorREF(ielemLOC,:) ; 
    % Coefficients evaluation polynomial 
    COEFFSpol = POLYINFO.COEFFSpolynomial{ielemLOC} ; 
    % 
    xLOC = (xLOC-coorREF)./LelemSCALING   ;  % Scaled coordinates 
    [Pevaluate,~]= CoordinateMatrixPolynomial(xLOC,VAR_SMOOTH_FE.ORDER_POLYNOMIALS)  ;
    Ninterpolation{ipoint} =  sparse(Pevaluate*COEFFSpol) ;  % Shape functions at the given points COORevaluate  
        
end

Ninterpolation = blkdiag(Ninterpolation{:}) ;  % Turn into a diagonal, sparse matrix... 

 
disp('*******************************************************************************')






if ~isempty(ELEMENTS_xNEW)
    disp('MESH CONTAINING THE NEW REDUCED SET OF POINTS');
    HYPERREDUCED_VARIABLES.WdomRED = wNEW ;
    HYPERREDUCED_VARIABLES.setElements = ELEMENTS_xNEW ;
  %  DATAIN.LABEL_NAME_PROJECT = DATA_GENGAUSS.NAMEFILE ;
%     DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'PLOT_REDUCED_GAUSS_IN_GID',0) ;
%     DATA_GENGAUSS = DefaultField(DATA_GENGAUSS,'SCALE_FACTOR_FOR_PRINTING_ECM_POINTS',0.1) ;
%     
%     if DATA_GENGAUSS.PLOT_REDUCED_GAUSS_IN_GID == 1
%         DATALOC.PRINT_MESH_ECM_POINTS_AS_POINTS = 1;
%         DATALOC.COORDINATES_POINTS_TO_PRINT = xNEW ;
%         DATALOC.COORDINATES_POINTS_TO_PRINT = xNEW ;
%         DATALOC.SCALE_FACTOR_FOR_PRINTING_ECM_POINTS = DATA_GENGAUSS.SCALE_FACTOR_FOR_PRINTING_ECM_POINTS ;
%     end
%     
%     
    DATALOC = DefaultField(DATALOC,'PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS',3) ; % =    1; % 0--> NONE % 1 -->Like time ... % 2 = Separated meshes
    
  %  if DATA_GENGAUSS.PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS == 0
   %     DATAIN = PrintingGid_ECMpoints(DATAIN,DATA_REFMESH,HROMVAR) ;
   % elseif DATA_GENGAUSS.PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS == 1
    %    DATAIN = PrintingGid_ECMhistory(DATAIN,DATA_REFMESH,HROMVAR,HISTORY_ITERATIONS) ;
    if DATALOC.PLOT_SELECTED_ELEMENTS_ALONG_ITERATIONS ==3
        PrintingGid_ECMpoints_LARGE(MESH,DATALOC,HISTORY_ITERATIONS) ;
    else
        error('Option not implemented')
    end
    
end

