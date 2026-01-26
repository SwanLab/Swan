function [ListElements,IndGaussLoc ]=  large2smallREP(IndGauss,ngaus)
%--------------------------------------------------------------------------
% function [ListElements, IndGaussLoc] = large2smallREP(IndGauss, ngaus)
%
% PURPOSE:
%   Given a list of global Gauss point indices (e.g., in a vectorized mesh),
%   this function maps them to:
%     1) Their corresponding *element number* (from which the Gauss point originates),
%     2) Their *local Gauss point index* within that element.
%
% DESCRIPTION:
%   - It assumes that each element has `ngaus` Gauss points.
%   - For each global Gauss point index in `IndGauss`, it computes:
%       - The associated element (via ceiling of index/ngaus),
%       - The local position within the element (1 to ngaus).
%
% INPUTS:
%   - IndGauss : [N x 1] Vector of global Gauss point indices.
%   - ngaus    : Scalar, number of Gauss points per element.
%
% OUTPUTS:
%   - ListElements : [K x 1] List of unique element indices involved (sorted by appearance).
%   - IndGaussLoc  : [N x 1] Local Gauss point index (1 to ngaus) within each corresponding element.
%
% EXAMPLE:
%   >> [ListElements, IndGaussLoc] = large2smallREP([1 2 3 4 9 20], 8)
%   ListElements = [1; 2; 3]
%   IndGaussLoc  = [1; 2; 3; 4; 1; 4]
%
% AUTHOR:
%   Joaquín A. Hernández, UPC-CIMNE
%--------------------------------------------------------------------------

if nargin == 0
    ngaus = 8 ; 
    IndGauss =[1 2 3 4 9 20 ] ; 
    
   
end

DIV = IndGauss/ngaus ;
ListElements = ceil(DIV) ; 
ListElements = ListElements(:) ;  

IndGaussINI = ListElements*ngaus - ngaus  ; 
IndGaussLoc = IndGauss(:)-IndGaussINI; 

ListElements = unique(ListElements,'stable') ; 
