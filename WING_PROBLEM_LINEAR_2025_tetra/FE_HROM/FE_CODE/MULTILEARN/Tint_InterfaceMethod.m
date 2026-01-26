function [Tcond,Tr ,cBC_d, cBC_r,BasisINT   ] = Tint_InterfaceMethod(BasisUdef,nDOM,BasisRdef,...
    DATAINM,BasisUrb,BasisRrb,alphaBC,NODESbound,ndim,COORref,uBAR,Hdr,Hrd,Hrr,...
    CNref,SingVal_disp,DATAIN,DATAOUT,Ui_Si,SingVal_reac,Ui_Si_reac)

if nargin == 0
    load('tmp2.mat')
end


DATA = [] ;

[BasisINT,f1,f2] = BasisInterfaceQ(BasisUdef,SingVal_disp,DATAIN,COORref,CNref,DATA,BasisRdef,DATAOUT,...
    Ui_Si,SingVal_reac,Ui_Si_reac) ;



%%% Component Td
BasisR = [ BasisRdef{1}] ;
ISTD = 1 ; 
Td =  Tmatrix(BasisR,f1,f2,nDOM,BasisINT,ISTD,DATAIN) ;
%%% Component Tr
BasisR = [ BasisRrb] ;
ISTD = 0 ;
Tr =  Tmatrix(BasisR,f1,f2,nDOM,BasisINT,ISTD,DATAIN) ;

Tcond = Td -Tr*(Hrr\Hrd) ;

%
%%% Boundary terms %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALPHA = zeros(1,2) ;
if alphaBC(1,1) ==1
    ALPHA(1) = 1;
end
if alphaBC(end,3) ==1
    ALPHA(2) = 1 ;
end
BasisR = [BasisRdef{1}];
cBC_d = zeros(size(BasisR,2)*nDOM,1) ;
if ALPHA(1) == 1
    COV= BasisR(f1,:)'*uBAR{1,1} ;
    IND = 1:size(COV,1) ;
    cBC_d(IND) = COV ;
end
if ALPHA(2) == 1
    COV= BasisR(f2,:)'*uBAR{end,3} ;
    IND = (size(cBC_d,1)-size(COV,1)+1):size(cBC_d,1) ;
    cBC_d(IND) = COV ;
end
BasisR = [BasisRrb];
cBC_r = zeros(size(BasisR,2)*nDOM,1) ;
if ALPHA(1) == 1
    COV= BasisR(f1,:)'*uBAR{1,1} ;
    IND = 1:size(COV,1) ;
    cBC_r(IND) = COV ;
end
if ALPHA(2) == 1
    COV= BasisR(f2,:)'*uBAR{end,3} ;
    IND = (size(cBC_r,1)-size(COV,1)+1):size(cBC_r,1) ;
    cBC_r(IND) = COV ;
end


%
% if ~isempty(Tint)
% nrows = size(Tint,1);
% nZ = nrows - size(Covf1,1) ;
% MATZEROS = sparse(nZ,size(Covf1,2)) ;
% else
%     MATZEROS = [] ;
% end


%
% NODESfaces = NODESbound.PLANE ;
% % f1 = small2large(NODESfaces{1},ndim) ;
% % f2 = small2large(NODESfaces{3},ndim) ;
%
% % Matrix BasisINTrb
% % -----------------------------------------
% f1 = NODESfaces{1} ; % Nodes involved face 1
% COOR_FACE = COORref(f1,:) ;
% COORrefPOINT = sum(COOR_FACE,1)/size(COOR_FACE,1); % Center of gravity
% COORrel = bsxfun(@minus,COOR_FACE',COORrefPOINT')'; % Relative coordinates
% BasisINTrb = ConstructBasisRigidBody(COORrel,DATAINM) ; %
% %%%%% BasisINTdef
% f1 = small2large(NODESfaces{1},ndim) ;
% f2 = small2large(NODESfaces{3},ndim) ;
% Xa = [BasisUdef{1}(f1,:),BasisUdef{1}(f2,:)] ;
% Xa = Xa -BasisINTrb*(BasisINTrb\Xa) ;
% [BasisINTdef S]= SVDT(Xa,0) ;
% % We take
% % % warning('Amend this ...')
% % % BasisINTdef = BasisINTdef(:,1:size(BasisUdef{1},2)) ;
%
% BasisINT = [BasisINTrb,BasisINTdef] ;

%
% Dbc = cell(2,1) ;
%
% if ALPHA(1) == 1
%     Dbc{1} = [Covf1;MATZEROS] ;
% else
%      Dbc{1} = [] ;
% end
% if ALPHA(2) == 1
%     Dbc{2} = [MATZEROS;Covf2] ;
% else
%      Dbc{2} = [] ;
% end
%
% %%% Boundary terms (reactions)
% BasisRintRB = BasisINTrb ;
% Xa = [BasisRdef{1}(f1,:),BasisRdef{1}(f2,:)] ;
% Xa = Xa -BasisRintRB*(BasisINTrb\Xa) ;
% [BasisRintDEF S]= SVDT(Xa,1e-6) ;
%
% BasisRint = [BasisRintRB,BasisRintDEF] ;
%
% JbcALL = -BasisINT'*BasisRint;
%
% %
% if ALPHA(1) == 1
%     cBC{1} = -BasisRint'*uBAR{1,1} ;
%     Jbc{1} = JbcALL ;
% else
%      cBC{1}= [] ;
%      Jbc{1} = [] ;
% end
% if ALPHA(2) == 1
%     cBC{2} = -BasisRint'*uBAR{end,3} ;
%       Jbc{2} = JbcALL ;
% else
%        cBC{2}= [] ;
%           Jbc{2} = [] ;
% end
% %%
% 
%
