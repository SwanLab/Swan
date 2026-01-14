function [G,uBAR,DOFr,DOFm,AREA,R,Fpnt] = PRESCRIBED_FLUCmixed_JOINTS_2faces_fun(F_A,F_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)

if nargin == 0
    load('tmp0.mat')
end

% DATA = DefaultField(DATA,'FLUCTUATION_MODES_BENDING',[]) ;
% Ub = DATA.FLUCTUATION_MODES_BENDING ;   % Not succesful...
load(DATA.FLUCTUATIONS_BOUNDARY_NameWS,'RIGID_BODY_MATRIX','FLUCTUATION_MATRIX','MASS_MATRIX_GEOMETRIC','FLUCTUATION_STIFF_COEFF')  ;

ndim = 3;
AREA = zeros(2,1) ;
R  = cell(2,1) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FACE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iface=1 ;
idomain = 1;
BasisINTfluct = FLUCTUATION_MATRIX{iface} ;  % Fluctuation modes (axial, torsion, .....)
RA = RIGID_BODY_MATRIX{iface} ;   % Rigid body modes
M = MASS_MATRIX_GEOMETRIC{iface} ;
TESTS = fieldnames(BasisINTfluct) ;

if iscell(RA)
    RA = RA{1} ; 
    M = M{1} ; 
end


U = [] ;
for itest = 1:length(TESTS)
    U  = [U, BasisINTfluct.(TESTS{itest})] ;
end
Util = U ;
MA = sparse(size(M,1)*ndim,size(M,2)*ndim) ;
for idim =1:ndim
    INDLOC =idim:ndim:(length(M)*ndim);
    MA(INDLOC,INDLOC) = M ;
end
nodesf = DOMAINVAR.NODES_faces12{idomain,iface} ;
[DOFA,rA,lA,JA] =   FlucMixed_Face(nodesf,Util,...
    RA,MA) ;
UbarA = MA*U ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FACE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iface=2 ;
idomain = size(DOMAINVAR.NODES_faces12,1);
BasisINTfluct = FLUCTUATION_MATRIX{iface} ;  % Fluctuation modes (axial, torsion, .....)
RB = RIGID_BODY_MATRIX{iface} ;   % Rigid body modes
M = MASS_MATRIX_GEOMETRIC{iface} ;
if iscell(RB)
    RB = RB{1} ; 
    M = M{1} ; 
end


TESTS = fieldnames(BasisINTfluct) ;
U = [] ;
for itest = 1:length(TESTS)
    U  = [U, BasisINTfluct.(TESTS{itest})] ;
end

Util = [RB,U]  ;
%nodesf = DOMAINVAR.NODES_faces12{idomain,iface} ;
MB = sparse(size(M,1)*ndim,size(M,2)*ndim) ;
for idim =1:ndim
    INDLOC =idim:ndim:(length(M)*ndim);
    MB(INDLOC,INDLOC) = M ;
end
UbarB = MB*U ;
nodesf = DOMAINVAR.NODES_faces12{idomain,iface} ;
[DOFB,rB,lB,JB] =   FlucMixed_Face(nodesf,Util,...
    RA,MB) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATRIX GAMMA
% This matrix depends on F_A and F_B, as well on the way in which  the
% fluctuations modes are sorted (typically: axial, torsion, bendy, bendz, shearz, sheary)
NameTestFE = DATA.NameTest  ;
kA = FLUCTUATION_STIFF_COEFF{1} ;
kB = FLUCTUATION_STIFF_COEFF{2} ;
GAMMA = zeros(1,length(TESTS))
for itest = 1:length(TESTS)
    nameCOEFF = TESTS{itest} ;
    GAMMA(itest) = kB.(nameCOEFF)/kA.(nameCOEFF) ;
    switch NameTestFE
        case 'sbend_y'
            switch nameCOEFF
                case 'shear_z'
                    L_m  = F_B(3)/F_B(5) ;
                    GAMMA(itest)  = GAMMA(itest)*(1+L_m)  ;
            end
        case 'sbend_z'
            switch nameCOEFF
                case 'shear_y'
                    L_m  = F_B(2)/F_B(6) ;
                    GAMMA(itest)  = GAMMA(itest)*(1+L_m)  ;
            end
            
    end
end
%% Recomputation of matrix UbarB  (including GAMMA )

UbarB = bsxfun(@times,UbarB',GAMMA')' ;

%%% CANONIC FORMAT: matrix G
%---------------------------
%
% METHOD_INCLUDING_GAMMA = 0 ;
%
% if METHOD_INCLUDING_GAMMA == 1
% Matrices L and H
%  \L \defeq \ident - \Ubar_{Al}^{-T}\Ubar_{Ar}^T   \J_A
L = UbarA(lA,:)'\(UbarA(rA,:)'*JA) ;
L = speye(size(L))-L ; %
%  \H  \defeq -\Ubar_{Al}^{-T}\Ubar_{Br}^T \J_B  +\Ubar_{Al}^{-T}\Ubar_{Bl}^T
H = -UbarA(lA,:)'\(UbarB(rB,:)'*JB)+UbarA(lA,:)'\(UbarB(lB,:)') ;
LinvH = L\H ;
DOFr = [DOFA(rA); DOFA(lA); DOFB(rB)] ;  % Slave DOFS
DOFm   = DOFB(lB) ; % Master DOFS
% --------------------------------------------------------
% MATRIX G
%G =  \coltres{-\J_A (\L^{-1}\H)}{\L^{-1}\H}{-\J_B}}}
G =  [-JA*LinvH;
    LinvH
    -JB] ;
% else
%     DOFr = [DOFA(rA); DOFB(rA) ] ;
%     DOFm = [DOFA(lA); DOFB(lB) ] ;
%     G = sparse(length(DOFr),length(DOFm));
%
%     inirow = 1;
%     finrow = size(JA,1) ;
%     inicol = 1;
%     fincol = size(JA,2);
%     G(inirow:finrow,inicol:fincol)  = -JA ;
%
%     inirow = finrow+1;
%     finrow = inirow+size(JB,1)-1 ;
%     inicol = fincol+1;
%     fincol = inicol+size(JB,2)-1;
%     G(inirow:finrow,inicol:fincol)  = -JB ;
%
%
% end

uBAR = zeros(length(DOFr),1) ;


%% EXTERNAL FORCES
% ----------------------------------
RbarB = MB*RB ;
Fpnt = zeros(size(COOR,1)*ndim,1) ;
b_B = (RbarB'*RB)\F_B ;
Fpnt(DOFB) = RbarB*b_B ;




