function [QrotNEW,norm_angle_rads] = UpdateRotationsEIFEM(OPERFE,Qrot,VAR,DATA,iterROTATIONS) ;
% Compute incremental rotations, and update rotation matrices (EIFEM, co-rot)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/02_PURE_ROTATION.mlx
% JAHO, 30-OCT-2024, Wednesday, 7:20, Balmes 185, Barcelona
% ---------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
% The amplitudes of the rigid body modes for all the EIF elements can be computed   as follows 
% \aRB  \DiagC{\PdownsRB} \Delta \dCqLOC
 aRB= OPERFE.D_PdownsRB*VAR.Delta_dCqLOC; 
 if DATA.MESH.ndim == 2
     aRBrot = aRB(3:3:end) ; 
     %  \thetaROT =  \DiagC{\lambdaLENall}^{-1}  \aRBrot
     thetaROT = aRBrot./OPERFE.lambdaLEN ; 
    % thetaROT = -aRBrot./OPERFE.lambdaLEN ; 
%      disp('Borrar estooooooo')
%      thetaROT = 1/180*pi*[1,1]';
     
     norm_angle_rads = norm(thetaROT)/length(thetaROT); 
     % Incremental rotation 
     IncQrot = zeros(size(Qrot)) ; 
%      ROWS = cell(1,2) ; 
%      COLS = cell(1,2) ; 
%      ROWS{1} = 1:2:size(IncQrot,1) ;
%      ROWS{2} = 2:2:size(IncQrot,1) ;
%      COLS{1}  = 1 ; 
%      COLS{2}  = 2 ; 
     ndim = 2; 
     IncQrot(1:ndim:size(IncQrot,1),1) = cos(thetaROT) ; 
     IncQrot(2:ndim:size(IncQrot,1),1) = sin(thetaROT) ; 
     IncQrot(1:ndim:size(IncQrot,1),2) = -sin(thetaROT) ; 
     IncQrot(2:ndim:size(IncQrot,1),2) = cos(thetaROT) ; 
     
     %\DiagC{\QrotALL } \leftarrow  \DiagC{\QrotALL }  \DiagC{\Delta \QrotALL(\thetaROT) } 
     % QrotNEW = Q*IncQrot ; 
%      
%      Q = Qrot(1:2,:) ; 
%      dQ = IncQrot(1:2,:) ; 
%      Qnew = Q*dQ; 


QrotNEW = MultiplyMatrixBlocks(Qrot,IncQrot) ; 
%      
%      QrotNEW = zeros(size(IncQrot)) ; 
%      nrows = size(QrotNEW,1) ; 
%      for iiiLOC = 1:ndim
%          iii = iiiLOC:ndim:nrows ;
%          for jjj = 1:ndim
%              for kkkLOC = 1:ndim
%                  kkk = kkkLOC:ndim:nrows;
%                  QrotNEW(iii,jjj) =   QrotNEW(iii,jjj)  + Qrot(iii,kkkLOC).*IncQrot(kkk,jjj) ;
%              end
%          end
%      end
     
     
     
 elseif DATA.MESH.ndim == 3
     aRBrot = reshape(aRB,6,[]) ; 
     aRBrot = aRBrot(4:6,:) ; 
     error('To be implemented')
 else
     error('Option not implemented')
 end 
 
 
 
formatSpec = '%10.5e'; 

 % disp(['iter=',num2str(iter),'; resf =',num2str(norm_res,formatSpec),' ; resf_ad=',num2str(criter_f,formatSpec),' res_disp_abs = ',num2str(normd,formatSpec)]) ;
  
  STRING_show = ['****** iterROT=%d  norm_angle =',formatSpec,'  strain_energy=',formatSpec ,'   \n'] ; 
  fprintf(1,STRING_show,[iterROTATIONS,norm_angle_rads,VAR.STRAIN_ENERGY]) ;