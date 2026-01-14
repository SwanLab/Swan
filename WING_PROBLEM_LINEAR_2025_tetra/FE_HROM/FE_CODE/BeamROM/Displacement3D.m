% function [DISP3D,DISP3D_lateral]= Displacement3D(DATAROM,MESH1D,qDEF,qRB,DATAIN,DATA_REFMESH) 
% 
% 
% if nargin ==0 
%     load('tmp1.mat')
% end
% 
% DISP3D = [] ; 
% DISP3D_lateral = [] ; 
% 
% 
% 
% DOMAINS_POSTPROCESS = DATAIN.DOMAINS_POSTPROCESS ; 
%  
% 
% qDEF = cell2mat(qDEF(DOMAINS_POSTPROCESS)') ; 
% qRB = cell2mat(qRB(DOMAINS_POSTPROCESS)') ; 
% 
% BasisUdef = DATAROM.BasisUdef ; 
% BasisUrb = DATA_REFMESH.BasisUrb; 
% 
% DISP3D = BasisUdef*qDEF + BasisUrb*qRB ; 
% 
% %% STRESSES 
% % ---------
%  
% 
% 
%  
% 
% 
% % 
% % V = DATAROM.BasisInt ;  % Interface modes matrix
% % ndim = size(V,2) ; % Number of entries for each node
% % nnode = size(MESH1D.COOR,1) ; % Number of 1D nodes (interfaces)
% % nelem = size(MESH1D.CN,1)  ; % Number of elements (slices)
% % nnodeE = 2; % Number of nodes per element (number of interfaces per element)
% % 
% % nmodesU = size(DATAROM.KdomRED,1) ; % Numbmer of displacement modes 
% % nmodesR = size(DATAROM.Hqr,2 ); % Number of reaction modes
% %  
% % F = cell2mat(fextDOMred) ;  % Reduced external froces
% % F = reshape(F,nmodesU,[]) ; 
% % 
% % rDEF = cell2mat(rDEF) ; 
% % rDEF = reshape(rDEF,nmodesR,[]) ;  % Amplitude self-equilibrated reactions
% % 
% % KHinv = DATAROM.KdomRED\DATAROM.Hqr ; 
% % qDEF = KHinv*rDEF  +DATAROM.KdomRED\F ; % Amplitude deformational displacements
% % 
% % %qDEFmatr = qDEF; 
% % %%%%%%%%%%%
% % qDEF = mat2cell(qDEF,size(qDEF,1),ones(1,size(qDEF,2))) ; 
% % qDEF = qDEF' ; 
% % 
% % %%% RIGID BODY AMPLITUDES 
% % f1 = DATAROM.f1 ;  % DOFs face 1
% % f2 = DATAROM.f2 ;  % DOFs face 2
% % f = [f1;f2] ; 
% % R_f = DATA_REFMESH.BasisUrb(f,:) ; % Rigid body modes faces 1 and 2 
% % 
% % T_1 = DATAROM.BasisInt'*DATAROM.BasisRdef(f1,:) ; 
% % T_2 = DATAROM.BasisInt'*DATAROM.BasisRdef(f2,:) ; 
% % 
% % RRinv = inv(R_f'*R_f) ; 
% % C_q = -RRinv*(R_f'*DATAROM.BasisUdef(f,:)) ; 
% % C_a{1} = RRinv*T_1 ; 
% % C_a{2} = RRinv*T_2 ; 
% %  
% % qRB = cell(nelem,1) ; 
% % for ielem = 1:nelem
% %     
% %     CNlocNOD = MESH1D.CN(ielem,:) ;
% %     % CNloc = Nod2DOF(CNlocNOD,ndim) ;
% %     qRB{ielem} = C_q*qDEF{ielem};
% %     for inode = 1:nnodeE
% %         NODE = CNlocNOD(inode) ;
% %         DOFS = Nod2DOF(NODE,ndim) ;
% %         qRB{ielem} = qRB{ielem} + C_a{inode}*a(DOFS)  ; 
% %     end
% %     
% % end
% % 
% % % qDEF = qDEFmatr; 
% % % qRB  =cell2mat(qRB); 
% % 
% % 
% % 
% % end
% % 
% %  