function   [xNEW,wNEW,POLYINFO,VAR_SMOOTH_FE,Ninterpolation,DATA_ECM,AUXVAR]...
    =  ConstrainedCECM...
    (W,COORg,HYPERREDUCED_VARIABLES,DATA_ECM,VAR_SMOOTH_FE)  ;

%(xFE,wFE,MESH,DATAIN,HYPERREDUCED_VARIABLES,DATA_ECM,VAR_SMOOTH_FE)  ;
if nargin == 0
    load('tmp1.mat')
end

VAR_SMOOTH_FE.COORg  = COORg ;
% Indexes integration points
z = HYPERREDUCED_VARIABLES.setPoints;
% REduced weights
w = HYPERREDUCED_VARIABLES.WdomRED ;
coorECM  = COORg(z,:) ;
% Exact Integral
DATA_ECM.ExactIntegral =  HYPERREDUCED_VARIABLES.PHI'*W ;

VAR_SMOOTH_FE.BasisIntegrand = HYPERREDUCED_VARIABLES.PHI ;
VAR_SMOOTH_FE.setPoints = z ;


% Approach after 2nd-January-2022
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/
% ContinuousEmpiricalCubatureM/Paper_hernandez2021ecm/07_BEAMS_NOGRAD/ControlledVariables.mlx

[xNEW,wNEW,DATA_ECM,VAR_SMOOTH_FE,POLYINFO,AUXVAR] = ConstrainedCECM_local(COORg,w,DATA_ECM,coorECM,VAR_SMOOTH_FE) ;

%  AUXVAR.PATHiniPOINTS = PATHiniPOINTS;
% AUXVAR.ELEMENTS_TO_PLOT = ELEMENTS_TO_PLOT;
% AUXVAR.HISTORY = HISTORY_ITERATIONS;


VAR_SMOOTH_FE.IndexECMini = z ;
%VAR_SMOOTH_FE.Bhred_interp =  BdomRED_interp ;

%ELEMENTS_xGAUSS= POLYINFO.ELEMENTS_CONTAINING_xNEW ;

VAR_SMOOTH_FE.ELEMENTS_xNEW =  POLYINFO.ELEMENTS_CONTAINING_xNEW ;
VAR_SMOOTH_FE.COORiniECM  = AUXVAR.HISTORY.POINTS{1}; % iNITIAL SET OF ECM POINTS
VAR_SMOOTH_FE.WEIGHTSiniECM  = AUXVAR.HISTORY.WEIGHTS{1}; % iNITIAL SET OF ECM POINTS
VAR_SMOOTH_FE.ELEMENTSiniECM  = AUXVAR.HISTORY.ELEMENTS_CONTAINING_POINTS{1}; % iNITIAL SET OF ECM POINTS

% In order to perform interpolations of any Gauss variable, it is necessary
% to have at one's dispossal the coefficients of the shape functions as
% well as the scaling factos of ELEMENTS_xNEW

%INTERPOLATION_INFO.COEFFSpolynomial = POLYINFO.COEFFSpolynomial(ELEMENTS_xNEW);
%INTERPOLATION_INFO.SCALING_VARIABLES.LENGTH = POLYINFO.SCALING_VARIABLES.LENGTH(ELEMENTS_xNEW);
%INTERPOLATION_INFO.SCALING_VARIABLES.LENGTH = POLYINFO.SCALING_VARIABLES.LENGTH(ELEMENTS_xNEW);
%INTERPOLATION_INFO.SCALING_VARIABLES.coorREF = POLYINFO.SCALING_VARIABLES.coorREF(ELEMENTS_xNEW,:);


% -------------------------------------
% INTERPOLATION FUNCTIONS
% --------------------------------------
if DATA_ECM.Method_Evaluation_Basis_Integrand == 1
    Ninterpolation = cell(length(POLYINFO.ELEMENTS_CONTAINING_xNEW ),1) ; %
    for ipoint = 1:length(POLYINFO.ELEMENTS_CONTAINING_xNEW )
        % Loop over elements
        % Coordinate of the point under consideration
        xLOC = xNEW(ipoint,:) ;
        ielemLOC = POLYINFO.ELEMENTS_CONTAINING_xNEW (ipoint) ;
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
else
    Ninterpolation = [] ;
end


disp('*******************************************************************************')





