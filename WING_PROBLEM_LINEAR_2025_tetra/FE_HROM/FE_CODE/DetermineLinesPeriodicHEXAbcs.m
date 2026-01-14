function [NODESln PLANES]= DetermineLinesPeriodicHEXAbcs(NODESpl)

%  See HEXAG_2.jpg

%

if nargin == 0
    load('tmp.mat')
end


%if length(NODESpl) == 6 
NODESln = cell(8,1) ; 
PLANES = cell(8,1) ; 

PLANES{1} = [1 5] ;
PLANES{2} = [2 5] ;
PLANES{3} = [3 5] ; 
PLANES{4} = [4 5] ;

PLANES{5} = [1 6] ;
PLANES{6} = [2 6] ;
PLANES{7} = [3 6] ; 
PLANES{8} = [4 6] ;

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




for i = 1:length(PLANES) ; 
    NODESln{i} = intersect(NODESpl{PLANES{i}(1)}, NODESpl{PLANES{i}(2)}) ; 
end