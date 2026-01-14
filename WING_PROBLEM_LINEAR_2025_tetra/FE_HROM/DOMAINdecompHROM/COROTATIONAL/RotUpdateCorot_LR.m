function  Qrot_new =   RotUpdateCorot_LR(OPERFE,Qrot,DATA,delta_dL,LboolCallQ,Delta_dR,kiter,VAR,DOFl,DOFr,PdownsRBcoupROTc_i)
% Update of incremental rotation matrix
% JAHO, 15-Feb-2025, Saturday, 16:19.  Balmes 185, Barcelona
% ---------------------------------------------------------------
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/05_COROT_SSLR_LSSR.mlx
if nargin == 0
    load('tmp1.mat')
end

Delta_dC = zeros(size(VAR.DISP)) ;
Delta_dC(DOFl) = delta_dL ;

% FIRST ITERATION
if kiter ==1
    Delta_dC(DOFr) = Delta_dR ;  % THIS DOES NOT SEEN TO HAVE A GREAT EFFECT 
end



% Local Coordinates
% Delta \dCloc = \DiagC{\QrotALL}^T  \Delta \dC
Delta_dCloc_allELEM = LboolCallQ*Delta_dC ;


if DATA.MESH.ndim == 2
    %\Delta \thetaROTloc = \DiagC{\PdownsRBlROT}  \Delta \dClocE{all}  + \dClocINCRE
    Delta_thetaROTloc = OPERFE.D_PdownsRBlROT*(Delta_dCloc_allELEM  + VAR.dClocINCRE) ;
    EXP_spin_Delta_thetaROT = zeros(size(Qrot)) ;
    EXP_spin_Delta_thetaROT(1:2:end,1) = cos(Delta_thetaROTloc) ;
    EXP_spin_Delta_thetaROT(2:2:end,1) = sin(Delta_thetaROTloc) ;
    EXP_spin_Delta_thetaROT(1:2:end,2) = -sin(Delta_thetaROTloc) ;
    EXP_spin_Delta_thetaROT(2:2:end,2) = cos(Delta_thetaROTloc) ;
    
    
    %
    %    dQrot_new = MultiplyMatrixBlocks(dQrot,EXP_spin_Delta_thetaROT) ;
    
else
     % 3D PROBLEMS 
     ndim = 3; 
     Delta_thetaROT = zeros(size(OPERFE.D_PdownsRBlROT,1),1) ; 
     
     for idim =  1:ndim 
         Delta_thetaROT(idim:ndim:end) = PdownsRBcoupROTc_i{idim}*(Delta_dCloc_allELEM  + VAR.dClocINCRE) ;
     end
     
     EXP_spin_Delta_thetaROT  =RodriguesFormula(Delta_thetaROT) ; 
     
% %      % Method 2  ---Gives the same.
%      Delta_thetaROTloc = OPERFE.D_PdownsRBlROT*(Delta_dCloc_allELEM  + VAR.dClocINCRE) ;
%      Delta_thetaROT_new = MultiplyMatrixVectorBlocks(Qrot,Delta_thetaROTloc) ;
    
end

    Qrot_new = MultiplyMatrixBlocks(EXP_spin_Delta_thetaROT,Qrot) ;
