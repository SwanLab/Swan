function  [z_mst, z_slv,w_mst,w_slv,DATA_regress_eta_der] = DiscreteECM_adaptWEIGHTSfst(A,z,w,DATA_interp,wFE)
%--------------------------------------------------------------------------
% DiscreteECM_adaptWEIGHTSfst
%
% PURPOSE:
%   Improved and accelerated version of `DiscreteECM_adaptWEIGHTS` for
%   computing the nonlinear master/slave mapping in the Empirical Cubature
%   Method (ECM) when the integrand lies on a low-dimensional nonlinear
%   manifold (typically 1D).
%
%   This version preserves the theoretical foundations and methodology
%   fully described in the original `DiscreteECM_adaptWEIGHTS` (see
%   references below), but introduces computational optimizations:
%     - Avoids repeated anonymous function calls when evaluating η(q),
%       η′(q), and η″(q).
%     - Builds and stores a precompiled evaluation structure
%       (`DATA_regress_eta_der`) via `BsplinesLeastSquares_fastECM` for
%       fast repeated evaluations in the online stage.
%
% METHOD OVERVIEW (as in the original version):
%   1. Compute exact integral from full data (A'*wFE) and ECM approximation.
%   2. Select a single "master" ECM point — the most explanatory in terms
%      of representing the others.
%   3. Treat all other ECM points as "slaves".
%   4. Apply weighted SVD to the slave data.
%   5. Fit a spline-based nonlinear mapping from master to slave values.
%   6. Use η(q) to reconstruct slave contributions from the master point.
%
% INPUTS:
%   - A           : Cell array or matrix of integrand snapshots 
%                   (M × Nsnap, M = Gauss points).
%   - z           : Indices of selected ECM points (length p).
%   - w           : Weights for the ECM points (length p).
%   - DATA_interp : Struct with interpolation settings:
%                     • LocalToleranceECMnonlinear_masterPOINTS
%                     • NSAMPLES, spline order, INCLUDE_SECOND_DERIVATIVES, etc.
%   - wFE         : Full finite element weights (M × 1), used for exact 
%                   reference integral computation.
%
% OUTPUTS:
%   - z_mst             : Index of the chosen master ECM point.
%   - z_slv             : Indices of slave ECM points.
%   - w_mst             : Weight of master ECM point.
%   - w_slv             : Weights of slave ECM points.
%   - DATA_regress_eta_der : Struct containing the precompiled spline-based 
%                            mapping and its derivatives for η(q).
%
% REFERENCES (original theory and method description):
%   - J.A. Hernández (2025), “Digression: on ECM applied to nonlinear manifolds”.
%   - Bravo et al. (2024), SAW-ECM, *International Journal for Numerical Methods
%     in Engineering*, 125(24), e7590.
%
% REMARKS:
%   - The adaptive ECM is designed for integrands on smooth nonlinear
%     manifolds where a single point can represent others via a nonlinear
%     function.
%   - This FAST variant is a drop-in replacement for the original in the
%     offline stage, producing identical results but faster online evaluation.
%
% DATE & PLACE OF MODIFICATION:
%   11-Aug-2025, Molinos Marfagones, Cartagena.
%
% AUTHOR:
%   Joaquín A. Hernández Ortega (UPC/CIMNE)
%   Comments updated by ChatGPT, 12-Aug-2025
%--------------------------------------------------------------------------
 
%--------------------------------------------------------------------------

if nargin == 0
    load('tmp1.mat')
    % DATA_interp.METHOD_SELECT_MASTER_ECM_POINTS = 'MAXIMUM_WEIGHT' ;
     
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

 

DATA_regress_eta_der = BsplinesLeastSquares_fastECM(DATA_interp, A_mst,  V_s', U_s, S_s);
A_slv = feval(DATA_regress_eta_der.nameFunctionEvaluate,A_mst,DATA_regress_eta_der) ; 
Int_AdaptiveECM =A_mst'*w_mst +  A_slv*w_slv ;


ERRO_ECMnew =  norm(Int_AdaptiveECM-ExactInt(indMST))/norm(ExactInt) ;
disp(['Error adaptive ECM = ',num2str(ERRO_ECMnew*100),' %'])


if  norm(ERRO_ECMstand-ERRO_ECMnew)/norm(ERRO_ECMstand) > 100.0 
    
    warning(['Either the fitting is not accurate, or simply there are no functions relating master/slave values (check the fitting)',  char(10),...
        '  Launch again the analysis  and select two master points'])
end

 