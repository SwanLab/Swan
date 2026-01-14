function [G,uBAR,DOFr,DOFm,AREA,R] = ...
    PRESCRIBED_FLUCTUATIONS_JOINTS2faces_COARSE(a_A,a_B,DOMAINVAR,COOR,CONNECTb,TypeElementB,...
    DATA,NAMECOARSE)

if nargin == 0
    load('tmp.mat')
end

load(DATA.FLUCTUATIONS_BOUNDARY_NameWS,'RIGID_BODY_MATRIX','FLUCTUATION_MATRIX','MASS_MATRIX_GEOMETRIC','CENTROIDS')  ;

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
[DOFrA,DOFmA,GA,uA,R{iface}] ...
    =  FLuctuBC_onefaceCOARSE(iface,DOMAINVAR,COOR,CONNECTb,TypeElementB,BasisINTfluct,...
    RIGID_BODY_AMPLITUDE,BasisINTrb,M,NAMECOARSE,CENTROIDS{iface}) ;
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

[DOFrB,DOFmB,GB,uB,R{iface}] ...
    =  FLuctuBC_onefaceCOARSE(iface,DOMAINVAR,COOR,CONNECTb,TypeElementB,BasisINTfluct,...
    RIGID_BODY_AMPLITUDE,BasisINTrb,M,NAMECOARSE,CENTROIDS{iface}) ;

DOFr = [DOFrA; DOFrB] ;
DOFm = [DOFmA; DOFmB] ;


G = sparse(length(DOFr),length(DOFm)) ;
uBAR = zeros(length(DOFr),1) ;

iini = 1;
ifin = size(GA,1) ;
iiniC = 1;
ifinC = size(GA,2) ;

G(iini:ifin,iiniC:ifinC) = GA ;
uBAR(iini:ifin) = uA  ;

iini = ifin+1;
ifin = iini+size(GB,1)-1 ;
iiniC = ifinC+1;
ifinC = iiniC+size(GB,2)-1 ;

G(iini:ifin,iiniC:ifinC) = GB ;
uBAR(iini:ifin) = uB ;


