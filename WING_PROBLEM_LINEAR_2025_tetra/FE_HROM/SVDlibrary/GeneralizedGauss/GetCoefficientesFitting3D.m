function DATAFITTING = GetCoefficientesFitting3D(xINI,G,DATALOC)
% POLYNOMIAL FITTING OF THE SVD DECOMPOSITION OF EACH ROWS OF G
% See GetCoefficientesFitting3D_aux.mlx
% JAHO, 20-Nov-2021
if nargin == 0
    load('tmp.mat')
    DATALOC.PLOT_functions_SVD_decompositionFITTING = 1;
end

DATALOC=  DefaultField(DATALOC,'TOL_SVDfitting',1e-6) ;  % Tolerance for the SVD (z level)

DATALOC=  DefaultField(DATALOC,'TOL_SVDfitting_xy',0) ;  % Tolerance for the SVD (xy level)


DATALOC=  DefaultField(DATALOC,'PLOT_MODES_G',0) ;  % Plot the functions
DATALOC = DefaultField(DATALOC,'PLOT_functions_SVD_decompositionFITTING',0);

% Loop over columns of PHI  --- Transponse of G

NFUN = size(G,1) ;
% Cell arrays to store spline coefficients
DATAFITTING.SPLINE_COEFFS_G_X = cell(NFUN,1)  ;  % Y functions
DATAFITTING.SPLINE_COEFFS_G_Y = cell(NFUN,1)  ;  % Y functions
DATAFITTING.SPLINE_COEFFS_G_Z = cell(NFUN,1)  ;  % X functions

DATAFITTING.DER_SPLINE_COEFFS_G_X = cell(NFUN,1)  ;  % Y functions
DATAFITTING.DER_SPLINE_COEFFS_G_Y = cell(NFUN,1)  ;  % Y functions
DATAFITTING.DER_SPLINE_COEFFS_G_Z = cell(NFUN,1)  ;  % X functions


DATAFITTING.SPLINE_COEFFS_SingularValues_Z = cell(NFUN,1)  ; % Singular values
DATAFITTING.SPLINE_COEFFS_SingularValues_XY = cell(NFUN,1)  ; % Singular values


disp('Computing SVD fitting parameters ...')
x = DATALOC.XYgrid{1} ;
y = DATALOC.XYgrid{2} ;
z = DATALOC.XYgrid{3} ;

% Loop over number of function
for iMODE = 1:NFUN
    disp(['G(',num2str(iMODE),',:)']) ;
    Gi = G(iMODE,:) ;
    lx = length(x) ;
    ly = length(y) ;
    % Reshaping Gi so that each row of C has the
    % same z coordinate
    C = reshape(Gi',lx*ly,[])';
    % To prove the above, check:
    %zzINIc  = reshape(xINI(:,3),lx*ly,[])';  % (z coordinate  matrix)
    % ---------------------------------------
    % SVD
    epsilonZ = DATALOC.TOL_SVDfitting;
    DATASVD.RELATIVE_SVD = 1;
    [Wz,Sz,VbarXY] = RSVDT(C,epsilonZ,[],0,DATASVD) ;
    
    
    rz = length(Sz) ; % Number of z modes
    % Spline coefficients, and nested SVDs
    Wz_splineCOEFF = cell(rz,1) ;
    Wz_splineCOEFF_der = cell(rz,1) ;
    
    Ux_splineCOEFF = cell(rz,1) ;
    Ux_splineCOEFF_der = cell(rz,1) ;
    
    Vy_splineCOEFF = cell(rz,1) ;
    Vy_splineCOEFF_der = cell(rz,1) ;
    
    Sxy =  cell(rz,1) ;
    
    for j = 1:rz
        Wz_splineCOEFF{j} =  spline(z,Wz(:,j)) ;
        Wz_splineCOEFF_der{j} = ppDer(Wz_splineCOEFF{j}) ;
        TOLfitting = DATALOC.TOL_SVDfitting_xy ;
        % Apply the same functionused in GetCoefficeintsFitting2D.m
        [Sxy{j},Ux_splineCOEFF{j},Ux_splineCOEFF_der{j},Vy_splineCOEFF{j},Vy_splineCOEFF_der{j}] =...
            SplineCoefficientsGetXY(x,y,VbarXY(:,j),TOLfitting) ;
    end
    
    DATAFITTING.SPLINE_COEFFS_G_Z{iMODE} = Wz_splineCOEFF  ;
    DATAFITTING.DER_SPLINE_COEFFS_G_Z{iMODE} = Wz_splineCOEFF_der  ;
    
    DATAFITTING.SPLINE_COEFFS_G_X{iMODE} = Ux_splineCOEFF  ;
    DATAFITTING.DER_SPLINE_COEFFS_G_X{iMODE} = Ux_splineCOEFF_der  ;

    DATAFITTING.SPLINE_COEFFS_G_Y{iMODE} = Vy_splineCOEFF  ;
    DATAFITTING.DER_SPLINE_COEFFS_G_Y{iMODE} = Vy_splineCOEFF_der  ;
    
    DATAFITTING.SPLINE_COEFFS_SingularValues_Z{iMODE} = Sz; 
    DATAFITTING.SPLINE_COEFFS_SingularValues_XY{iMODE} = Sxy; 
    
%     if DATALOC.PLOT_functions_SVD_decompositionFITTING == 1
%         
%         PlotSVDmodesFitting(DATALOC,iMODE,S,Ux,Vy,Ux_splineCOEFF,Vy_splineCOEFF,x,y) ;
%     end
    
    
end
