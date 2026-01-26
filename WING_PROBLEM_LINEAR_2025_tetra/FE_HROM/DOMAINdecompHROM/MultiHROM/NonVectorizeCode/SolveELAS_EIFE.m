function [d,  stressGLO,TRANSF_COORD,stressesREF   ] = SolveELAS_EIFE(K,Fb,Ftrac,dR,DOFr,COOR,...
    CN,TypeElement,DATA,M,NameFileMesh,PROPMAT,MaterialType,TRANSF_COORD,Bmat,WEIGHTSinteg) ;
% SOLUTION OF THE EIFE equations
% JAHO, 12-MARCH-2023
% ----------------------
if nargin == 0
    load('tmp.mat')
end
nnode = size(COOR,1); ndim = size(COOR,2); nelem = size(CN,1); nnodeE = size(CN,2) ;     %
% Solution of the system of FE equation
% Right-hand side
%dbstop('20')
F = Fb + Ftrac ;
% Set of nodes at which temperature is unknown
DOFl = 1:nnode*ndim ;
DOFl(DOFr) = [] ;  

% dL =  K^{-1}*(Fl .Klr*dR)
%dbstop('27')
dL = K(DOFl,DOFl)\(F(DOFl)-K(DOFl,DOFr)*dR) ;
% Vector of   displacements
d = zeros(nnode*ndim,1) ;
d(DOFl)= dL ;
d(DOFr) = dR ;
% Reaction forces
React = zeros(size(d)) ;
%dbstop('33')
React(DOFr) = K(DOFr,:)*d -F(DOFr) ;

%

% %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MomentsResult =  ComputeResultantReactions(DOFr,React,COOR,ndim); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% COmputation of heat flux vector at each gauss point
 disp('Computation of stress and strains at each element (maximum)')
   [stressGLO,TRANSF_COORD,stressesREF]   = StressStrains_EIFE(COOR,CN,TypeElement,d,DATA,PROPMAT,MaterialType,TRANSF_COORD,Bmat,WEIGHTSinteg) ;
% strainGLO = [] ; 
% stressGLO = [] ; posgp = [] ; 
posgp = [] ; 
DATA = DefaultField(DATA,'NAME_FILE_STORE','DATAMK'); 

nameWS = [DATA.NAME_FILE_STORE,'_',DATA.NAME_INPUT_DATA,'.mat'] ;  
save(nameWS,'K','M','d','DOFl','DOFr','COOR','CN','TypeElement','posgp','NameFileMesh','Ftrac')
