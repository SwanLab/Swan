function [G,uBAR,DOFr,DOFm,AREA,R] = ...
    PRESCRIBED_FLUCTUATIONS_JOINTS2faces_fun(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA)

if nargin == 0
    load('tmp.mat')
end

% DATA = DefaultField(DATA,'FLUCTUATION_MODES_BENDING',[]) ;
% Ub = DATA.FLUCTUATION_MODES_BENDING ;   % Not succesful...
load(DATA.FLUCTUATIONS_BOUNDARY_NameWS,'RIGID_BODY_MATRIX','FLUCTUATION_MATRIX','MASS_MATRIX_GEOMETRIC')  ;

ndim = 3;
AREA = zeros(2,1) ; 
R  = cell(2,1) ; 
%%%% FACE 1
iface=1 ;
BasisINTfluct = FLUCTUATION_MATRIX{iface} ; 
BasisINTrb = RIGID_BODY_MATRIX{iface} ; 
RIGID_BODY_AMPLITUDE = a_A ; 
M = MASS_MATRIX_GEOMETRIC{iface} ; 
if iscell(M)
    M = M{1} ;
    BasisINTrb = BasisINTrb{1} ; 
end
[DOFrA,DOFmA,JA,uA,AREA(iface),R{iface},MstA] ...
    =  FLuctuBC_oneface(iface,DOMAINVAR,COOR,CONNECTb,TypeElementB,BasisINTfluct,...
    RIGID_BODY_AMPLITUDE,BasisINTrb,M) ; 
%%%% FACE 2
iface=2 ;
BasisINTfluct = FLUCTUATION_MATRIX{iface} ; 
BasisINTrb = RIGID_BODY_MATRIX{iface} ; 
RIGID_BODY_AMPLITUDE = a_B ;
M = MASS_MATRIX_GEOMETRIC{iface} ; 
if iscell(M)
    M = M{1} ;
    BasisINTrb = BasisINTrb{1} ; 
end

[DOFrB,DOFmB,JB,uB,AREA(iface),R{iface},MstB] ...
    =  FLuctuBC_oneface(iface,DOMAINVAR,COOR,CONNECTb,TypeElementB,BasisINTfluct,...
    RIGID_BODY_AMPLITUDE,BasisINTrb,M) ; 

DOFr = [DOFrA; DOFrB] ; 
DOFm = [DOFmA; DOFmB] ; 
  
 
G = sparse(length(DOFr),length(DOFm)) ;
uBAR = zeros(length(DOFr),1) ;

iini = 1;
ifin = size(JA,1) ;
iiniC = 1;
ifinC = size(JA,2) ;

G(iini:ifin,iiniC:ifinC) = JA ;
uBAR(iini:ifin) = uA  ;

iini = ifin+1;
ifin = iini+size(JB,1)-1 ;
iiniC = ifinC+1;
ifinC = iiniC+size(JB,2)-1 ;

G(iini:ifin,iiniC:ifinC) = JB ;
uBAR(iini:ifin) = uB ;


