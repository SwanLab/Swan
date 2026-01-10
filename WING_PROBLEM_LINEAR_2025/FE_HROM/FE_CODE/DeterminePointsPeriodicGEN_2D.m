function [NODESpnt LINES]= DeterminePointsPeriodicGEN_2D(NODESln)

 
%

if nargin == 0
    load('tmp.mat')
end

LINES = cell(4,1);
%if length(NODESpl) == 6 
NODESpnt = cell(4,1) ; 
 
LINES{1} = [1 2] ;
LINES{2} = [2 3] ;
LINES{3} = [1 4] ; 
LINES{4} = [3 4] ;
 


for i = 1:length(LINES) ; 
    NODESpnt{i} = intersect(NODESln{LINES{i}(1)}, NODESln{LINES{i}(2)}) ; 
end