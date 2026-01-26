function [NODESln PLANES]= DetermineLinesPeriodic(NODESpl)


%

if nargin == 0
    load('tmp.mat')
end


NODESln = cell(18,1) ; 
PLANES = cell(18,1) ; 

PLANES{1} = [1 7] ;
PLANES{2} = [2 7] ;
PLANES{3} = [3 7] ; 
PLANES{4} = [4 7] ;
PLANES{5} = [5 7] ;
PLANES{6} = [6 7] ; 
 
PLANES{7} = [1 8] ;
PLANES{8} = [2 8] ;
PLANES{9} = [3 8] ; 
PLANES{10} = [4 8] ;
PLANES{11} = [5 8] ;
PLANES{12} = [6 8] ;

PLANES{13} = [1 2] ;
PLANES{14} = [2 3] ; 
PLANES{15} = [3 4] ; 
PLANES{16} = [4 5] ; 
PLANES{17} = [5 6] ; 
PLANES{18} = [6 1] ; 

for i = 1:length(PLANES) ; 
    NODESln{i} = intersect(NODESpl{PLANES{i}(1)}, NODESpl{PLANES{i}(2)}) ; 
end