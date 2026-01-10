function    invJ= inverseVECTORIZED(Jinp,ndim,detJinp)
% Given a matrix Jinp = [Jinp_1; Jinp_2 ... ; Jinp_nelem]  (Jinp_i is ndim x ndim),
% inverseTRANSvectorize returns a matrix consisting of the    
%   inverse of each block matrix
%  invJ = [inv(Jinp_1'); inv(Jinp_2') ... ;inv(Jinp_nelem') ]
% ndim = 2,3
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 21-Oct-2021
if nargin == 0
    J1 = [3 5 6; 7 8 7 ; 0 0 3] ;  d1 =inv(J1) ;
    J2 = 4*[3 5 6; 7 8 7 ; 0 0 3] ;  d2 =inv(J2) ;
    J3 = 6*[3 5 6; 7 8 7 ; 0 0 3] ;  d3 =inv(J3) ;
    
        ndim = 3 ;

    J1 = [3 5 ; 7 8  ] ;  d1 =inv(J1) ;
    J2 = 4*[3 5 ; 8 8] ;  d2 =inv(J2) ;
    J3 = 6*[3 5 ;  8 3 ] ;  d3 =inv(J3) ;
    
    Jinp = [J1; J2; J3]; invJ_check = [d1;d2;d3] ;
    detJinp = [] ;
    ndim = 2 ;
end

if isempty(detJinp)
     detJinp= determinantVECTORIZE(Jinp,ndim) ; 
end

%
J = cell(ndim,ndim) ;
for i=1:ndim
    iglo = i:ndim:size(Jinp,1) ;
    for j=1:ndim
        jglo = j ;
        J{i,j} = Jinp(iglo,jglo) ;
    end
end
% ------------------------------------------------
    invJ = zeros(size(detJinp)) ;

if ndim ==2
      % i = 1
    % -----
    i= 1; iglo = i:ndim:size(Jinp,1) ;
    %---------------------------------
    invJ(iglo,1) = J{2,2}./detJinp ;
    invJ(iglo,2) = - J{1,2}./detJinp ;
    % i = 2
    % -----
    i= 2; iglo = i:ndim:size(Jinp,1) ;
    %---------------------------------
    invJ(iglo,1) =  -J{2,1}./detJinp ;
    invJ(iglo,2) =   J{1,1}./detJinp ;
    
elseif ndim == 3
    % i = 1
    % -----
    i= 1; iglo = i:ndim:size(Jinp,1) ;
    %---------------------------------
    invJ(iglo,1) = (J{2,2}.*J{3,3} - J{3,2}.*J{2,3})./detJinp ;
    invJ(iglo,2) = (J{3,2}.*J{1,3} - J{1,2}.*J{3,3})./detJinp ;
    invJ(iglo,3) = (J{1,2}.*J{2,3} - J{2,2}.*J{1,3})./detJinp ;
    % i = 2
    % -----
    i= 2; iglo = i:ndim:size(Jinp,1) ;
    %---------------------------------
    invJ(iglo,1) = (J{3,1}.*J{2,3} - J{2,1}.*J{3,3} )./detJinp ;
    invJ(iglo,2) = (J{1,1}.*J{3,3} - J{3,1}.*J{1,3} )./detJinp ;
    invJ(iglo,3) = (J{2,1}.*J{1,3} - J{1,1}.*J{2,3} )./detJinp ;
    % i = 3
    % -----
    i= 3; iglo = i:ndim:size(Jinp,1) ;
    %---------------------------------
    invJ(iglo,1) = (J{2,1}.*J{3,2} - J{3,1}.*J{2,2} )./detJinp ;
    invJ(iglo,2) = ( J{3,1}.*J{1,2} - J{1,1}.*J{3,2} )./detJinp ;
    invJ(iglo,3) = (J{1,1}.*J{2,2} - J{2,1}.*J{1,2})./detJinp ;
end

 


