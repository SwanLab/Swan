function [PhiMaster_nonl,PhiSlave_nonl,qMASTER_nonl,...
    qSLAVE_nonl,...
    DATA_evaluateTAU_and_DER] = Bsplines_EP_1stmode(Phi_non,Kll,DATA_interp,SNAPdisp_plast_L)
% =========================================================================
% BSPLINES_EP_1STMODE — Plastic Manifold Mapping via First Plastic Mode
% =========================================================================
% PURPOSE
%   Define the nonlinear plastic manifold using the FIRST plastic SVD mode as
%   the MASTER latent coordinate. All remaining plastic modes are SLAVES and
%   are regressed as functions of the master via least-squares B-splines. The
%   mapping is returned in a compact SVD form for robustness and efficiency.
%
% WORKFLOW (high level)
%   1) Split plastic basis:
%        PhiMaster_nonl = Phi_non(:,1),  PhiSlave_nonl = Phi_non(:,2:end).
%   2) Compute Kll-weighted coordinates for plastic snapshots:
%        qMASTER_nonl = PhiMaster_nonl' * Kll * SNAPdisp_plast_L
%        qSLAVE_nonl  = PhiSlave_nonl'  * Kll * SNAPdisp_plast_L
%      – Remove duplicate master entries; normalize sign so training range
%        is non-negative when appropriate (tension/compression handling).
%   3) (Optional) Symmetry/exploitation & preprocessing:
%      – If DATA_interp.ExploitSymmetryPlasticLatentVariable = true, the master
%        abscissa can be |qMASTER| (single-sided training).
%      – A monotonicity/spacing transform may be applied to reduce clustering.
%   4) Sampling strategy for regression:
%      – If DATA_interp.SubSampling_qPLASTIC_master = true, subsample along
%        qMASTER to obtain near-uniform coverage and avoid overfitting near
%        yield onset; otherwise use all snapshots.
%   5) Compact SVD of slave coordinates:
%        [U,S,V] = SVDT(qSLAVE_nonl(_SAMP))
%      – Represent slaves in a low-rank basis before spline fitting.
%   6) Least-squares B-spline regression:
%      – Fit each component of V (or V^T) as a spline of qMASTER using
%        DATA_interp.NSAMPLES, order_Bslines, ratio_NSAMPLES_knots, etc.
%      – Assemble callable decoder and derivatives: τ(q), τ′(q), τ″(q).
%
% INPUTS
%   Phi_non          : [nDOF × r_pl] plastic basis (Kll-orthonormal columns).
%   Kll              : Constrained stiffness matrix (energy-metric).
%   DATA_interp      : Struct with interpolation/fit controls:
%                        • NSAMPLES, order_Bslines, ratio_NSAMPLES_knots
%                        • SubSampling_qPLASTIC_master (logical)
%                        • ExploitSymmetryPlasticLatentVariable (logical)
%   SNAPdisp_plast_L : [|DOFl| × n_plastic] plastic snapshots (projected to DOFl,
%                      already orthogonalized w.r.t. elastic mode elsewhere).
%
% OUTPUTS
%   PhiMaster_nonl        : First plastic mode (master).
%   PhiSlave_nonl         : Remaining plastic slave modes.
%   qMASTER_nonl          : Kll-weighted coordinates along master.
%   qSLAVE_nonl           : Kll-weighted coordinates along slaves.
%   DATA_evaluateTAU_and_DER : Struct with function handles to evaluate
%                              τ(q), τ′(q), τ″(q) at runtime.
%
% PRACTICAL NOTES / GOTCHAS
%   • MASTER must be the first plastic SVD mode to guarantee an invertible
%     encoder and consistent stress evaluation downstream.
%   • Subsampling the master abscissa improves spline conditioning and reduces
%     bias from snapshot clustering.
%   • Use symmetry exploitation only when training is single-sided in sign and
%     the constitutive response is symmetric in that regime.
%
% DEPENDENCIES
%   SVDT, BsplinesLeastSquares_PLAST, BsplinesSUBSAMPL,
%   transformOddPointsScaled, DefaultField.
%
% VERSION / AUTHORSHIP
%   • 18-AUG-2025 — J.A. Hernández, Cartagena. First version based on
%     elastoplastic manifold with first-mode master variable.
%   • 07-NOV-2025 — Comments refreshed; symmetry/sign normalization and
%     subsampling clarified. Barcelona.
%   Author: Joaquín Alberto Hernández Ortega (JAHO) — jhortega@cimne.upc.edu
%   (This header was generated automatically by ChatGPT on 07-NOV-2025.)
% =========================================================================

% =========================================================================


if nargin == 0
    load('tmp1.mat')
%       DATA_interp.NSAMPLES = 300 ;
%       DATA_interp.ratio_NSAMPLES_knots = 0.8;
%       DATA_interp.SubSampling_qPLASTIC_master = true;  
    close all
end

% NumberOfPoints_masterRANGE = size(SNAPdisp_plast_L,2) ;
%   
%
% DATA_interp = DefaultField(DATA_interp,'NumberOfPoints_masterRANGE',NumberOfPoints_masterRANGE) ;
%



% This is the master mode (first mode)
iMASTER  = 1;
PhiMaster_nonl =    Phi_non(:,iMASTER) ;
iSLAVE = setdiff(1:size(Phi_non,2),iMASTER) ;
PhiSlave_nonl =  Phi_non(:,iSLAVE) ; % Slave modes, plastic
%
qSLAVE_nonl = PhiSlave_nonl'*(Kll*SNAPdisp_plast_L) ; % Slave coefficients
qMASTER_nonl = PhiMaster_nonl'*(Kll*SNAPdisp_plast_L) ; % Master coefficients

[qMASTER_nonl,aaa,bbb ]= unique(qMASTER_nonl) ; % REmove repeated indices

signTRAIN = sign(qMASTER_nonl) ; 
signTRAIN = unique(signTRAIN) ; 

if signTRAIN == -1
    qMASTER_nonl = -qMASTER_nonl ; 
    PhiMaster_nonl = -PhiMaster_nonl ; 
    signTRAIN = 1;  
end

DATA_interp.signTRAIN = signTRAIN; 

NumberSigns = length(signTRAIN) ; 
DATA_interp = DefaultField(DATA_interp,'ExploitSymmetryPlasticLatentVariable',false) ; 

if NumberSigns  == 2 && DATA_interp.ExploitSymmetryPlasticLatentVariable
    error('Train either in tension or compression ') 
elseif  DATA_interp.ExploitSymmetryPlasticLatentVariable
    qMASTER_nonl_input = abs(qMASTER_nonl) ; 
else
    qMASTER_nonl_input = qMASTER_nonl ; 
end


 

figure(10)
hold on
xlabel('Number of snapshot')
ylabel('qMASTER nonl')
plot(qMASTER_nonl,'Marker','*')
grid on


 y_transformed = transformOddPointsScaled(qMASTER_nonl, 1/2) ; 
 
 plot(y_transformed,'Marker','*')

 

qSLAVE_nonl = qSLAVE_nonl(:,aaa) ;

DATA_interp = DefaultField(DATA_interp,'SubSampling_qPLASTIC_master',true) ;

if ~DATA_interp.SubSampling_qPLASTIC_master
    
    [UU,SS,VV] = SVDT(qSLAVE_nonl) ;
    [DATA_evaluateTAU_and_DER] = BsplinesLeastSquares_PLAST(DATA_interp, qMASTER_nonl_input, VV', UU, SS) ;
    
else
  [DATA_evaluateTAU_and_DER] =  BsplinesSUBSAMPL(qMASTER_nonl_input,DATA_interp,qSLAVE_nonl) ; 
    
    
    
end


