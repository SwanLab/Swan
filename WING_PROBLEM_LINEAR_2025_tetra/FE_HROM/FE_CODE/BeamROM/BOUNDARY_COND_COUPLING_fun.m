function [G,uBAR,DOFr,DOFm,AREA,R,Fpnt] = ...
    BOUNDARY_COND_MIXED_fun(DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA,b_A_input,b_B_input)

if nargin == 0
    load('tmp.mat')
end

DATA = DefaultField(DATA,'NameWS_bending_displacements',[]) ;

if ~isempty(DATA.NameWS_bending_displacements)
    load(DATA.NameWS_bending_displacements,'dBENDING')  ;
else
    dBENDING =[] ;
end




ndim = 3;

INFOBENDING.dBENDING = dBENDING ;
if b_B_input(5) ~= 0
     INFOBENDING.ORDER_MODE = 1;
    INFOBENDING.b = [3,5] ;
elseif b_B_input(6) ~= 0
    % Shear test in the z direction
         INFOBENDING.ORDER_MODE = 2;
    INFOBENDING.b = [2,6] ;
end


%%%% FACE 1
% ----------

if isempty(dBENDING)
    [DOFA,AREA,R,M,Ub,QA,QB,G,uBAR,DOFr,DOFm,DOFB]...
        =  MIXED_CONDITIONcoup_face(DOMAINVAR,COOR,CONNECTb,TypeElementB,...
        INFOBENDING) ;
    
     
    
else
    [DOFA,AREA,R,M,Ub,QA,QB,G,uBAR,DOFr,DOFm,DOFB]...
        =  SIMPLE_BENDING_coup_face(DOMAINVAR,COOR,CONNECTb,TypeElementB,dBENDING,...
        DATA,INFOBENDING,b_B_input) ;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------------
%%%%% NODAL FORCES
% ----------------------------------
Fpnt = zeros(size(COOR,1)*ndim,1) ;
b_A = (R'*M*R)\b_A_input ;
b_B = (R'*M*R)\b_B_input ;
% Nodal forces FACE A
Fpnt(DOFA) = M*R*b_A ;
% Nodal forces FACE B
Fpnt(DOFB) = M*R*b_B ;




%
% %%%
%
% G = sparse(length(DOFr),length(DOFm)) ;

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
%
%
% %
% % METHOD = 1;
% %
% %
% % if METHOD ==1
% %     % Include constraint that Ub'[d_A ;d_B] = 0
% %       ZEROS_M = zeros(size(R,2),length(DOFB)) ;
% %
% %     A = [QA -QB
% %          R'*M ZEROS_M] ;
% %
% %
% %
% %     factor = 0.5 ;
% %
% %     A = [(Ub') factor*(Ub')
% %           (Ib(s,:))  (-Ib(s,:))
% %          (M*R)' ZEROS_M] ;
% %
% % %      A = [    (Ib(s,:))  (-Ib(s,:))
% % %          (M*R)' ZEROS_M] ;
% %
% %      b = [zeros(size(Ub,2),1);
% %            R(s,:)*(a_A-a_B)
% %           (R'*M*R)*a_A] ;
% %
% % %
% % %         b = [          R(s,:)*(a_A-a_B)
% % %           (R'*M*R)*a_A] ;
% %
% %       % Choose a set of linearly independent equations
% % %       A_aug = [A,b] ;
% % %      [~,INDEX] =  licols(A_aug') ;
% % %      A = A(INDEX,:) ;
% % %      b = b(INDEX) ;
% % %
% %      % Now select length(INDEX) linearly ind. columns from A
% %       [~,DOFr] =  licols(A) ;
% %       DOFr = DOFr(:) ;
% %
% %       DOFm = setdiff(1:size(A,2),DOFr) ;
% %       DOFm = DOFm(:) ;
% %       % Ar*dR + Am*dM = b --> dR = inv(Ar)*(b-Am*dM)
% %       G = -A(:,DOFr)\(A(:,DOFm)) ;
% %       G = sparse(G) ;
% %       uBAR = A(:,DOFr)\b ;
% %
% %       DOFAB = [DOFA; DOFB] ;
% %       DOFr = DOFAB(DOFr) ;
% %       DOFm = DOFAB(DOFm) ;
% %
% %
%
%
% %else
%
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
%
%
%
