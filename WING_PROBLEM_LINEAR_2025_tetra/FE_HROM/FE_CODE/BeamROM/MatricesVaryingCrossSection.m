function [A,D,a0] = MatricesVaryingCrossSection(GEOref,CROSS_SECTION,s,L) ;
%
%  Transformation matrices (scaling, varying cross-sections), see
%  BeamROM.pdf
%
%  x = a0 + A*X +X_1*D*X
%  where
%  a0  = [0 a2 a3]' ;
%  A = [e1  0  0 ;
%       b2  c2  0
%       b3   0  c3]
%  D = [0  0   0
%       0  d2  0
%       0   0    d3]
% JAHO, 25 Januyary -2018
%b  --------------------------

if nargin == 0
    load('tmp1.mat')
end


% Variations in x direction (scaling)
l =  GEOref.dX_0 ;
ds = abs(s(2)-s(1)) ;


e1 = ds/l;


a1 = 0 ;

% Variations in y directions
FIELDS = {'YMAX','YMIN'} ;
[a2,b2,c2,d2] = ConstantMatricesTransformation(s,l,L,GEOref,CROSS_SECTION,FIELDS) ;
% Variations in z directions

if isfield(GEOref,'ZMAX')
    FIELDS = {'ZMAX','ZMIN'} ;
    [a3,b3,c3,d3] = ConstantMatricesTransformation(s,l,L,GEOref,CROSS_SECTION,FIELDS) ;
    
    a0 = [a1;a2;a3] ;
    A = [e1 0 0
        b2 c2 0
        b3  0   c3] ;
    D = [0  0  0
        0  d2  0
        0  0  d3] ;
    
else
    
      a0 = [a1;a2] ;
    A = [e1 0 
        b2 c2 
        ] ;
    D = [0  0  
        0  d2   
        ] ;
    
    
end

