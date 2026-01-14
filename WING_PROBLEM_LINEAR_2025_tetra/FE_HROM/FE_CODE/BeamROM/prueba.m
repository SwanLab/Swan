clc
clear all
load('tmp,mat')
%  
% load(DATA.NameWS_bending_displacements,'dBENDING')  ;

dA = d(DOFA) ; 


% 
% ndim = 3;
% %%%% FACE 1
% % ----------
% iface=1 ;
% [DOFA,AREA,R,M,L_sp,hAB,J_rl,J_rp,b_As,r,l,s,p]...
%     =  FLuctBC_onefaceBENDING(iface,DOMAINVAR,COOR,CONNECTb,TypeElementB,dBENDING,...
%     a_A,a_B,DATA) ; 
% %%%% FACE 2
% % ----------
% iface=2 ;
% nodesfB = DOMAINVAR.NODES_faces12{end,iface} ;
% DOFB = small2large(nodesfB,ndim) ;   % DOFS face 1
%%%
% Slave and master DOFs
% DOFr = [DOFA(s(r)); DOFB(s(r)) ; DOFB(s(l))] ; 
% DOFm = [DOFA(s(l)); DOFA(p) ; DOFB(p)] ; 
% G = sparse(length(DOFr),length(DOFm)) ;
% uBAR = zeros(length(DOFr),1) ;
% 
% %-----------------------------------------------------
% % First row %-----------------------------------------
% % ----------------------------------------------------
% iini = 1;
% ifin = size(J_rl,1) ;
% % --- 1st column
% iiniC = 1;
% ifinC = size(J_rl,2) ;
% G(iini:ifin,iiniC:ifinC) = -J_rl ;
% % % --- 2nd column
% iiniC = ifinC + 1;
% ifinC = iiniC + size(J_rp,2)-1  ;
% G(iini:ifin,iiniC:ifinC) = -J_rp ;
% % Ind. term 
% uBAR(iini:ifin) = b_As  ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %-----------------------------------------------------
% % Second row %-----------------------------------------
% % ----------------------------------------------------
% iini = ifin+1;
% ifin = iini + size(J_rl,1)-1 ;
% % --- 1st column
% iiniC = 1;
% ifinC = size(J_rl,2) ;
% G(iini:ifin,iiniC:ifinC) = -J_rl ;
% % % --- 2nd column
% iiniC = ifinC + 1;
% ifinC = iiniC + size(J_rp,2)-1  ;
% G(iini:ifin,iiniC:ifinC) = -J_rp + L_sp(r,:) ;
% % % --- 3rd column
% iiniC = ifinC + 1;
% ifinC = iiniC + size(J_rp,2)-1  ;
% G(iini:ifin,iiniC:ifinC) = -L_sp(r,:)  ;
% % Ind. term 
% uBAR(iini:ifin) = b_As -hAB(r)  ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %-----------------------------------------------------
% % Third row %-----------------------------------------
% % ----------------------------------------------------
% iini = ifin+1;
% ifin = iini + length(l)-1 ;
% % --- 1st column
% iiniC = 1;
% ifinC = size(J_rl,2) ;
% G(iini:ifin,iiniC:ifinC) = eye(length(l),length(l)) ;
% % % --- 2nd column
% iiniC = ifinC + 1;
% ifinC = iiniC + length(p)-1  ;
% G(iini:ifin,iiniC:ifinC) =  L_sp(l,:) ;
% % % --- 3rd column
% iiniC = ifinC + 1;
% ifinC = iiniC + length(p)-1  ;
% G(iini:ifin,iiniC:ifinC) = -L_sp(l,:)  ;
% % Ind. term 
% uBAR(iini:ifin) =  -hAB(l)  ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% 
% 
% 
% 
%  
% 
% 
