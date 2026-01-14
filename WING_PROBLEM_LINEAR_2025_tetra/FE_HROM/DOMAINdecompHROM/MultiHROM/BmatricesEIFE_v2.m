function [Bmat,w,x,MATPRO,RECONS_PK2stresses,index_elements,BmatRED,Bgrad,PhiDEFelem]  = ...
    BmatricesEIFE_v2(PhiDEF,PdownsDEF,CECM_intforces,BstFE,DATAoffline,DATA,MATPRO,BasisStwo,MESH)
% B-matrix empirical interscale FE
% JAHO, 23-Feb-2023/1-MArch-2023
% Updated, improved version of BmatricesEIFE_v1.m
% See comments in /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/04_DistortedElement.mlx
% -------------------------------
if nargin == 0
    load('tmp.mat') % CECM
    %   DATAoffline.UseDECMpoints = 0;
end

% Location/weights integration points
% ----------------------------
%DATAoffline = DefaultField(DATAoffline,'UseDECMpoints',0) ;

if DATAoffline.UseDECMpoints.INTERNAL_FORCES == 1
   % DECM points 
    [Bmat,w,x,MATPRO,RECONS_PK2stresses,index_elements,BmatRED,Bgrad,PhiDEFelem]  = ... ; 
    BmatricesEIFE_v2_DECM(PhiDEF,PdownsDEF,CECM_intforces,BstFE,DATAoffline,DATA,MATPRO,BasisStwo,MESH) ; 
    
else
    % CECM points 
    [Bmat,w,x,MATPRO,RECONS_PK2stresses,index_elements,BmatRED,Bgrad,PhiDEFelem]  = ... ; 
    BmatricesEIFE_v2_CECM(PhiDEF,PdownsDEF,CECM_intforces,BstFE,DATAoffline,DATA,MATPRO,BasisStwo,MESH) ; 
    
end
