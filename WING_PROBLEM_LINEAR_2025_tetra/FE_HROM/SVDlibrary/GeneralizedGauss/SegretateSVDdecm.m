function  [Lambda,S,V,ETIME,ERRORsvd,RankMatrix,DATAOUT] =    SegretateSVDdecm(Xf, DATA )
if nargin == 0
    load('tmp.mat')
end
DATA = DefaultField(DATA,'TOLsvdXf_ElementsHighAccuracy',0) ; 

COLS_2 = DATA.PointsToBeIntegratedWithHighAccuracy; 

npoints = size(Xf,1) ; 
COLS_1 = setdiff(1:npoints,COLS_2) ; 
 

XfT_tol = cell(1,2) ; 
DATAsvd.EPSILON_GLO  =zeros(1,2) ; 

idom = 1;  % Standard domain 
XfT_tol{idom} = Xf(COLS_1,:)'; 
DATAsvd.EPSILON_GLO(1) = DATA.TOLsvdXf ; 

idom = 2;  % High-accuracy domain 
XfT_tol{idom} = Xf(COLS_2,:)'; 
DATAsvd.EPSILON_GLO(2) = DATA.TOLsvdXf_ElementsHighAccuracy ; 


[U,S,V,ETIME,ERRORsvd,RankMatrix,DATAOUT] =    RSVDqp(XfT_tol, DATAsvd.EPSILON_GLO ) ;


[~,INV_COLS] = sort([COLS_1,COLS_2']) ; 
Lambda = V(INV_COLS,:) ; 
V = U; 