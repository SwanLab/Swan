function [NODESpnt ]= DeterminePointPeriodicHEXA(NODESln)

% See Implementation.pdf

if nargin == 0
    load('tmp.mat')
end


%if length(NODESpl) == 6 
LINES = cell(8,1) ; 

LINES{1} = [1 5] ;
LINES{2} = [2 5] ;
LINES{3} = [3 5] ; 
LINES{4} = [4 5] ;

LINES{5} = [1 6] ;
LINES{6} = [2 6] ;
LINES{7} = [3 6] ; 
LINES{8} = [4 6] ;

% 
% else
%     
% NODESln = cell(8,1) ; 
% LINES = cell(8,1) ; 
% 
% LINES{1} = [1 5] ;
% LINES{2} = [2 5] ;
% LINES{3} = [3 5] ; 
% LINES{4} = [4 5] ;
% 
% LINES{5} = [1 6] ;
% LINES{6} = [2 6] ;
% LINES{7} = [3 6] ; 
% LINES{8} = [4 6] ;
%     
% end




for i = 1:length(LINES) ; 
    NODESln{i} = intersect(NODESpl{LINES{i}(1)}, NODESpl{LINES{i}(2)}) ; 
end