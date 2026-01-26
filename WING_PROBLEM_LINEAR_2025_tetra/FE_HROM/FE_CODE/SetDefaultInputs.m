function [DATA DOFm Gb MaterialType]= SetDefaultInputs(DATA)
% Set default values


DATAdef.PLOT.REACTIONS = 1 ;
DATAdef.CALCULATE_averageSTRESS = 0 ;
DOFm = [] ;
Gb = [] ;
DATAdef.VECTcode  = 1 ;   % Vectorize code
DATAdef.BCBformulation = 1 ;  % Formulation B^T C B
MaterialType = [] ;
DATAdef.NOVOIDS = 1 ;
DATAdef.TYPESOLVER = 0 ; % Options: "default": 0 / conjugated gradient (1)
DATAdef.niterCONJG = 1000 ; % Number of iterations  conjugated gradient
DATAdef.tolCONJG = 1e-6 ; % Tolerance
DATAdef.RECALCULATE_STIFFNESS = 1 ;
DATAdef.plotFLUCT  =0 ; % Plot fluctuations
DATAdef.STORE_STIFFNESS = 1 ; %Store stiffness matrix
DATAdef.REFERENCE_POINT = 'CENTER'; % Options: 'CENTER','CORNER' (FOR COMPUTING TOTAL DISPLACEMENTS )
DATAdef.SPECIAL_BC_SHEAR_ON = 0 ; % For laminates. Special study of shear stiffness factors 
DATAdef.ORDER_SHEAR_THEORY = 1 ;  % 1 or 3
DATAdef.CALCULATE_STRENGTH_LAM= 1 ; 
DATAdef.CALCULATE_AVGDENS = 0 ; 
DATAdef.CALCULATE_MASSMATRIX = 0 ; 


%   DATA.REFERENCE_POINT = 'CENTER';  % Options: 'CENTER','CORNER' (FOR COMPUTING TOTAL DISPLACEMENTS )
%     DATA.STORE_STIFFNESS = 0 ; %Store stiffness matrix
%     DATA.niterCONJG = 10000 ; % Number of iterations  conjugated gradient
%     DATA.tolCONJG = 1e-10 ; % Tolerance

fff =fieldnames(DATAdef) ;

for i=1:length(fff)
    fieldLOC =fff{i} ;
    if ~isfield(DATA,fieldLOC)
        DATA.(fieldLOC) = DATAdef.(fieldLOC) ;
    end
end