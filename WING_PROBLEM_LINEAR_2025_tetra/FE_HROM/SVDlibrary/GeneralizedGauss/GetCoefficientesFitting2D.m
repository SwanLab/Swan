function DATAFITTING = GetCoefficientesFitting2D(xINI,G,DATALOC)
% POLYNOMIAL FITTING OF THE SVD DECOMPOSITION OF EACH ROWS OF G
% JAHO, 20-Nov-2021
if nargin == 0
    load('tmp1.mat')
    DATALOC.PLOT_functions_SVD_decompositionFITTING = 1;
end

DATALOC=  DefaultField(DATALOC,'TOL_SVDfitting',1e-6) ;  % Tolerance for the SVD of each column
DATALOC=  DefaultField(DATALOC,'PLOT_MODES_G',0) ;  % Plot the functions
DATALOC = DefaultField(DATALOC,'PLOT_functions_SVD_decompositionFITTING',0);

% Loop over columns of PHI  --- Transponse of G

NFUN = size(G,1) ;
DATAFITTING.SPLINE_COEFFS_G_dim{1} = cell(NFUN,1)  ;
DATAFITTING.SPLINE_COEFFS_G_dim{2} = cell(NFUN,1)  ;
DATAFITTING.DER_SPLINE_COEFFS_G_dim{1} = cell(NFUN,1)  ;
DATAFITTING.DER_SPLINE_COEFFS_G_dim{2} = cell(NFUN,1)  ;
DATAFITTING.SPLINE_COEFFS_SingularValues = cell(NFUN,1)  ;

% Sall_numb = zeros(size(G,1),1) ;
% Errorx_all = Sall_numb ;
% Errory_all = Sall_numb ;
disp('Computing SVD fitting parameters ...')
x = DATALOC.XYgrid{1} ;
y = DATALOC.XYgrid{2} ;

% Loop over number of function
for iMODE = 1:NFUN
    disp(['G(',num2str(iMODE),',:)']) ;
    Gi = G(iMODE,:) ;
    TOL = DATALOC.TOL_SVDfitting;
    [S,Ux_splineCOEFF,Ux_splineCOEFF_der,Vy_splineCOEFF,Vy_splineCOEFF_der,Ux,Vy] =...
        SplineCoefficientsGetXY(x,y,Gi,TOL) ;
    
    
    if DATALOC.PLOT_functions_SVD_decompositionFITTING == 1
        
        PlotSVDmodesFitting(DATALOC,iMODE,S,Ux,Vy,Ux_splineCOEFF,Vy_splineCOEFF,x,y) ;
    end
    
    DATAFITTING.SPLINE_COEFFS_G_dim{1}{iMODE} = Ux_splineCOEFF  ;
    DATAFITTING.SPLINE_COEFFS_G_dim{2}{iMODE} = Vy_splineCOEFF ;
    DATAFITTING.DER_SPLINE_COEFFS_G_dim{1}{iMODE} = Ux_splineCOEFF_der  ;
    DATAFITTING.DER_SPLINE_COEFFS_G_dim{2}{iMODE} = Vy_splineCOEFF_der ;
    DATAFITTING.SPLINE_COEFFS_SingularValues{iMODE} = S ;
end

%
%
%
%
% %     DATALOC.NFUN = iMODE;
% %     DATALOC.PLOT_FUNCTIONS = DATALOC.PLOT_MODES_PHI ;
% %
%
%     [Uypoly,Vxpoly,Uypoly_der,Vxpoly_der,S,Errorx,Errory] = ...
%         SplineSVDinterpolation_poly(xx,yy,Fmat,DATALOC) ;
%
%     Sall_numb(iMODE) = length(S) ;
%
%     DATAFITTING.Uypoly{iMODE} = Uypoly;
%     DATAFITTING.Vxpoly{iMODE} = Vxpoly;
%     DATAFITTING.Uypoly_der{iMODE} = Uypoly_der;
%     DATAFITTING.Vxpoly_der{iMODE} = Vxpoly_der ;
%     DATAFITTING.S{iMODE} = S ;
%
%
%     if ~isempty(Errorx)
%         Errorx_all(iMODE) = Errorx;
%         Errory_all(iMODE) = Errory;
%         disp(['Fitting ERROR x =',num2str(Errorx)])
%         disp(['Fitting ERROR x =',num2str(Errory)])
%     end
% end
%
%
% figure(34)
% subplot(2,1,1)
% hold on
% xlabel('Number of PHI column')
% ylabel('Number of Singular Values ')
% title('POLYNOMIAL FITTING SVD')
% bar(Sall_numb)
%
% subplot(2,1,2)
% hold on
% xlabel('Number of PHI column')
% ylabel('Fitting error')
% plot(Errorx_all)
% plot(Errory_all)
%
%
%
%
% %
% % [DATAFITTING] = SVDinterpolation2D_FITTINGspline(DATA,G',DATA.XYgrid{1},DATA.XYgrid{2}) ;
% %
% % % % Reshaping the matrix of initial points
% % % % --------------------------------------
% % % x = xINI(:,1) ; y  = xINI(:,2) ;
% % % % Gi = sum_j lambda_j  U_j(x)  V_j(y)^T
% % % % Variations in x should occur along columns, and variations in y along rows
% % % MPOINTS = DATA.MPOINTS ;
% % % xGRID = reshape(x,MPOINTS(2),[])' ;
% % % yGRID = reshape(y,MPOINTS(2),[])' ;
% % %
% % %
% % % for imodes = 1:nmodes
% % %     % Loop over number of modes
% % %     % Reshape Gi
% % %     C = reshape(G(imodes,:),MPOINTS(2),[])' ;
% % % end
% % %
% % %
% % %
% % %
% % %
% % % DATAFITTING.spline_G = cell(1,nmodes) ; % Coefficients splines for each column
% % % DATAFITTING.spline_derG = cell(1,nmodes) ;          % Coefficients derivatives splines for each column
% % %  for i = 1:nmodes
% % %     DATAFITTING.spline_G{i} = spline(xINI,G(i,:)) ;
% % %     DATAFITTING.spline_derG{i} =  fnder( DATAFITTING.spline_G{i}  ) ;
% % % end