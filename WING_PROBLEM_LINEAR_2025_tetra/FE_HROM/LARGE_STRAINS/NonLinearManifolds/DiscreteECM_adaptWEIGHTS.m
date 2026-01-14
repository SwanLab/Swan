function  [z_mst, z_slv,w_mst,w_slv,etaNON,etaNONder,etaNONder2] = DiscreteECM_adaptWEIGHTS(A,z,w,DATA_interp,wFE)
%--------------------------------------------------------------------------
% function [z_mst, z_slv, w_mst, w_slv, etaNON, etaNONder, etaNONder2] = ...
%     DiscreteECM_adaptWEIGHTS(A, z, w, DATA_interp, wFE)
%
% PURPOSE:
%   Implements an adaptive variant of the Empirical Cubature Method (ECM)
%   for integrands that lie on a low-dimensional nonlinear manifold—
%   typically a 1D curve in high-dimensional space. The method reduces the
%   cost of numerical integration by reconstructing all but one of the ECM
%   points (called slaves) from a single “master” point via nonlinear interpolation.
%
%   The key idea is to:
%     - Identify the master ECM point (most representative or explanatory),
%     - Use a weighted SVD on the remaining slave points,
%     - Interpolate the slave contributions as a nonlinear function of the master value.
%
%   This approach is useful for reduced-order modeling, where integrands
%   change smoothly with respect to a few reduced coordinates.
%
% INPUTS:
%   A           : Cell array or matrix of integrand snapshots (size M × Nsnap),
%                 where M is the number of Gauss points.
%   z           : Vector of selected ECM indices (length p).
%   w           : ECM weights associated with indices in z (length p).
%   DATA_interp : Struct containing interpolation settings and tolerances:
%                   * LocalToleranceECMnonlinear_masterPOINTS – SVD tolerance
%                   * NSAMPLES, INCLUDE_SECOND_DERIVATIVES, etc. (for B-spline)
%   wFE         : Vector of full finite element weights (size M×1), used to
%                 compute the exact reference integral.
%
% OUTPUTS:
%   z_mst       : Index of the selected master ECM point (scalar).
%   z_slv       : Vector of indices for the slave ECM points (length p-1).
%   w_mst       : Scalar weight of the master point.
%   w_slv       : Vector of weights for the slave points.
%   etaNON      : Nonlinear reconstruction map handle:
%                 etaNON(u_mst) ≈ [u_slv_1; u_slv_2; ...]
%   etaNONder   : Derivative of etaNON (optional use).
%   etaNONder2  : Second derivative of etaNON (if enabled in DATA_interp).
%
% METHOD OVERVIEW:
%   1. Compute exact integral from full data (A'*wFE).
%   2. Compute ECM approximation with current points (A(z,:)' * w).
%   3. Select most explanatory point (master) using error-based criterion.
%   4. Partition remaining points as slaves.
%   5. Apply weighted SVD to slave data (A_slv) with slave weights (w_slv).
%   6. Fit spline-based interpolation: A_mst → projection coefficients.
%   7. Define nonlinear map etaNON to reconstruct slave values from master.
%   8. Estimate new ECM integral using etaNON and compare with the exact one.
%
% REFERENCES:
%   - J.A. Hernández (2025), “Digression: on ECM applied to nonlinear manifolds”
%   - Bravo et al. (2024), SAW-ECM, *International Journal for Numerical Methods in Engineering*, 125(24), e7590.
%
% EXAMPLE:
%   [zm, zs, wm, ws, eta, deta] = DiscreteECM_adaptWEIGHTS(A, z, w, interp_data, wFE);
%   u_mst     = A(zm,:);           % Master value at new sample
%   u_slv_hat = eta(u_mst);        % Predicted slave values
%
% REMARKS:
%   - ECM points are assumed to lie on a smooth nonlinear manifold.
%   - The first basis in the SVD of A_slv is typically sufficient due to low intrinsic dimension.
%   - The approximation is validated by reassembling the integral and comparing with full FE.
%
% AUTHOR:
%   Joaquín A. Hernández, UPC/CIMNE, 4 July 2025
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
    % DATA_interp.METHOD_SELECT_MASTER_ECM_POINTS = 'MAXIMUM_WEIGHT' ;
    DATA_interp.LocalToleranceECMnonlinear_masterPOINTS = 1e-2;
end

A = cell2mat(A) ;
ExactInt = A'*wFE;
A_z = A(z,:);
Int_ECM_stand = A(z,:)'*w ;

ERRO_ECMstand =  norm(Int_ECM_stand-ExactInt)/norm(ExactInt) ;
disp(['Error Standard ECM = ',num2str(ERRO_ECMstand*100),' %'])
 

[UU,SS,ZZ] = SVDT(A_z); 
if length(SS) ~= length(z)
    error('A_z must be full rank')
end


% NOW WE CONFRONT THE MOST DELICATE PHASE OF THE ENTIRE "ADAPTIVE" ECM: WE
% HAVE TO ANSWER THE QUESTION:  IS THERE ANY SINGLE POINT (THE MASTER
% POINT) such that we can express the remaining variables as nonlinear
% functions of the integrand at this master point ?

  
%DATA_interp = DefaultField(DATA_interp,'METHOD_SELECT_MASTER_ECM_POINTS','MOST_EXPLANATORY') ;

% switch DATA_interp.METHOD_SELECT_MASTER_ECM_POINTS
%     case 'MOST_EXPLANATORY'
total_error = BestPointExplanatory(A_z') ;
disp('Indexes points sorted by degree of "representativeness (explanatory capabilities)" ')
[AAA,indEXPLAN] = sort(total_error) ;
disp(num2str(indEXPLAN)) ;

iCHOSEN  = 1;
best_col_index = indEXPLAN(iCHOSEN) ;
%     case 'MAXIMUM_WEIGHT'
%         [AAA,best_col_index] = max(w) ;
% end



compl_Set = setdiff(1:size(A_z,1),best_col_index) ;

A_mst = A_z(best_col_index,:) ;
A_slv = A_z(compl_Set,:) ;
z_mst = z(best_col_index) ;
w_mst = w(best_col_index) ;
z_slv = z(compl_Set);
w_slv = w(compl_Set) ;

PLOT_EVOLUTION_Az = 0 ;

if  PLOT_EVOLUTION_Az == 1
    disp('WARNING: The following graphical representation is useless (or misleading) in more than one latent dimension')
    figure(3094)
    hold on
    xlabel('Snapshot ')
    ylabel('Internal work density')
    numberECMpoints = size(A_z,1);
    % numberECMpoints = min(numberECMpoints, size(A_z,1)) ;
    title(['Internal work density first ',num2str(numberECMpoints),' ECM points' ])
    
    for ipoints = 1:numberECMpoints
        y = A_z(ipoints,:)/norm(A_z(ipoints,:)) ;
        plot(y,'DisplayName',num2str(ipoints))  ;
    end
    
    legend show
    
    
end

% 5-08-2025:  Why did we introduce this in the first place?
% See discussion in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/112_NonLIN_ROM_RBF/05_meta_AECM.mlx

% = DefaultField(DATA_interp,'METHOD_SORT_ENTRIES_MASTER_POINTS','monotonic') ;

% switch DATA_interp.METHOD_SORT_ENTRIES_MASTER_POINTS
%     case 'monotonic'
[A_mst,indMST] = sort(A_mst) ;
A_slv = A_slv(:,indMST) ;
%     case 'none'
%
%
%
% end

%











DATA_interp = DefaultField(DATA_interp,'LocalToleranceECMnonlinear_masterPOINTS',1e-5);


TOL = DATA_interp.LocalToleranceECMnonlinear_masterPOINTS;
M  =diag(w_slv) ;
DATA_SVD.RELATIVE_SVD = 1;
DATA_SVD.TOL = TOL ;
% This is the weighted SVD,
[U_s,S_s,V_s] = WSVDT(A_slv,M,DATA_SVD) ;


PLOT_BEFORE_SPLINES = 1;

if  PLOT_BEFORE_SPLINES == 1
    figure(463)
    hold on
    title('Right singular vectors SVD Amaster as function of Aslave')
    xlabel('A_{mst}')
    ylabel('V_{slv}')
    for iii  =1:size(V_s,2)
        plot(A_mst,V_s(:,iii))
    end
end


DATA_interp.Nfigures_base = 450 ;

DATA_interp = DefaultField(DATA_interp,'Legend_var_x','A_{mst}') ;
DATA_interp = DefaultField(DATA_interp,'Legend_var_y','\eta  (A_{slv})') ;

[etaNON, etaNONder,etaNONder2] = BsplinesLeastSquares_ECM(DATA_interp, A_mst, V_s', U_s, S_s) ;
%[etaNON, etaNONder,etaNONder2] = BsplinesLeastSquares(DATA_interp, A_mst, V_s', U_s, S_s) ;

%[etaNON, etaNONder,nREDcoor] = SplineInterpLOC_ECM(DATA_interp, A_mst, V_s', U_s, S_s) ;


%Int_AdaptiveECM =  etaNON(A_mst)'*[w_mst; w_slv] ;
Int_AdaptiveECM =A_mst'*w_mst +   etaNON(A_mst)'*w_slv ;


ERRO_ECMnew =  norm(Int_AdaptiveECM-ExactInt(indMST))/norm(ExactInt) ;
disp(['Error adaptive ECM = ',num2str(ERRO_ECMnew*100),' %'])


%if  norm(ERRO_ECMstand-ERRO_ECMnew)/norm(ERRO_ECMstand) > 10.0 
    
    warning(['Either the fitting is not accurate, or simply there are no functions relating master/slave values (check the fitting)',  char(10),...
        '  Launch again the analysis  and select two master points'])
%end

 