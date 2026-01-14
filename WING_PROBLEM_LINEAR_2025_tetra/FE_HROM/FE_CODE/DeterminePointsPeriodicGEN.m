function [NODESpnt LINES]= DeterminePointsPeriodicGEN(NODESln)

%  See HEXAG_2.jpg

%

if nargin == 0
    load('tmp.mat')
end

LINES = cell(8,1);
%if length(NODESpl) == 6 
NODESpnt = cell(8,1) ; 
 
LINES{1} = [1 10] ;
LINES{2} = [1 9] ;
LINES{3} = [9 5] ; 
LINES{4} = [10 5] ;

LINES{5} = [4 3] ;
LINES{6} = [3 11] ;
LINES{7} = [11 7] ; 
LINES{8} = [7 12] ;

 


% 
% else
%     
% NODESln = cell(8,1) ; 
% PLANES = cell(8,1) ; 
% 
% PLANES{1} = [1 5] ;
% PLANES{2} = [2 5] ;
% PLANES{3} = [3 5] ; 
% PLANES{4} = [4 5] ;
% 
% PLANES{5} = [1 6] ;
% PLANES{6} = [2 6] ;
% PLANES{7} = [3 6] ; 
% PLANES{8} = [4 6] ;
%     
% end




for i = 1:length(LINES) ; 
    NODESpnt{i} = intersect(NODESln{LINES{i}(1)}, NODESln{LINES{i}(2)}) ; 
end