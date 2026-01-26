function [Bmat,w,x,MATPRO,RECONS_PK2stresses,index_elements,BmatRED,Bgrad,PhiDEFelem]  = ...
    BmatricesEIFE_BUB(PhiDEF,PdownsDEF,CECM_intforces,BstFE,DATAoffline,DATA,MATPRO,BasisStwo,MESH)
% Modification of BmatricesEIFE_v2.m 
% Incorporation of bubble modes
% JAHO, 30-Sept-2023, Satudary, Barcelona (Farga by Secrets)
% -------------------------------
if nargin == 0
    load('tmp.mat') % CECM
    %   DATAoffline.UseDECMpoints = 0;
end

% Location/weights integration points
% ----------------------------
%DATAoffline = DefaultField(DATAoffline,'UseDECMpoints',0) ;

if DATAoffline.UseDECMpoints.INTERNAL_FORCES == 1
   % DECM points , no continuous version of the ECM ---just select among
   % Gauss points
    [Bmat,w,x,MATPRO,RECONS_PK2stresses,index_elements,BmatRED,Bgrad,PhiDEFelem]  = ... ; 
    BmatricesEIFE_v2_DECM(PhiDEF,PdownsDEF,CECM_intforces,BstFE,DATAoffline,DATA,MATPRO,BasisStwo,MESH) ; 
    
else
    % CECM points 
    [Bmat,w,x,MATPRO,RECONS_PK2stresses,index_elements,BmatRED,Bgrad,PhiDEFelem]  = ... ; 
    BmatricesEIFE_v2_CECM(PhiDEF,PdownsDEF,CECM_intforces,BstFE,DATAoffline,DATA,MATPRO,BasisStwo,MESH) ; 
    
end
