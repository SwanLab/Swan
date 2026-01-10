function dQrot_new =   RotUpdateCorot2(OPERFE,dQrot,DATA,ResCOMPall,Delta_dC,D_QrotALL) 
% Update of incremental rotation matrix
% Latex notation:: /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTexAPPV.tex 
%/home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/04_UPDrotUNC.mlx
% JAHO, 3-Dec-2024, Balmes 185, Barcelona 
% ---------------------------------------------------------------

% \Delta \aRBrot =   -\DiagC{\mFICTq(\cgreen{\dQrot})^{-1}} \Par{  \ResCOMPall   +  \DiagC{\YcmpD} \DiagC{\QrotALL}^T      \cblue{\Delta \dC}}
 if nargin == 0
     load('tmp1.mat')
 end 

if DATA.MESH.ndim == 2
    
    % \mFICTq = \DiagC{\Zfict} \dQrotLOC(:)  
    % dQrotLOC = dQrot
    dQrotLOC_vector = zeros(size(dQrot,1)*size(dQrot,2),1) ; 
    dQrotLOC_vector(1:4:end) = dQrot(1:2:end,1) ; 
    dQrotLOC_vector(2:4:end) = dQrot(2:2:end,1) ;
    dQrotLOC_vector(3:4:end) = dQrot(1:2:end,2) ;
    dQrotLOC_vector(4:4:end) = dQrot(2:2:end,2) ;
    mFICTq = OPERFE.D_Zfict*dQrotLOC_vector ; 
    Delta_aRBrot = ResCOMPall + OPERFE.D_YcmpD*(D_QrotALL'*(OPERFE.LboolCall*Delta_dC)) ; 
    Delta_aRBrot = -Delta_aRBrot./mFICTq ; 
    % \Delta\thetaROT =  \DiagC{\lambdaLENall}^{-1}   \Delta \aRBrot
    Delta_thetaROT = Delta_aRBrot./OPERFE.lambdaLEN ; 
    EXP_spin_Delta_thetaROT = zeros(size(dQrot)) ; 
    
    EXP_spin_Delta_thetaROT(1:2:end,1) = cos(Delta_thetaROT) ;
    EXP_spin_Delta_thetaROT(2:2:end,1) = sin(Delta_thetaROT) ;
    EXP_spin_Delta_thetaROT(1:2:end,2) = -sin(Delta_thetaROT) ;
    EXP_spin_Delta_thetaROT(2:2:end,2) = cos(Delta_thetaROT) ;
    
   
   % \dQrot \leftarrow  \DiagC{\dQrot}   \EXP{\spin{\Delta \thetaROT}}
   
   dQrot_new = MultiplyMatrixBlocks(dQrot,EXP_spin_Delta_thetaROT) ; 
%    
%    dQrot_new  =  zeros(size(dQrot))  ;   
%    ndim =2 ; 
%    for irowsLOC = 1:ndim 
%        irows = irowsLOC:ndim:size(dQrot_new,1);
%        for jcols = 1:ndim 
%            for kcolsLOC  = 1:ndim
%                 kcols = kcolsLOC:ndim:size(dQrot_new,1);
%            dQrot_new(irows,jcols) = dQrot_new(irows,jcols)  + dQrot(irows,kcolsLOC).*EXP_spin_Delta_thetaROT(kcols,jcols) ; 
%            end
%        end
%    end
     
else 
    
    error('Option not implemented')
    
end