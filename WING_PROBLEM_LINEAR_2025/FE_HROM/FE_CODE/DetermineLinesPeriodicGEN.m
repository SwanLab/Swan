function [NODESln PLANES]= DetermineLinesPeriodicGEN(NODESpl)

%  See HEXAG_2.jpg

%

if nargin == 0
    load('tmp.mat')
end

PLANES = cell(12,1) ; 
%if length(NODESpl) == 6 
NODESln = cell(12,1) ; 
 
PLANES{1} = [1 5] ;
PLANES{2} = [2 5] ;
PLANES{3} = [3 5] ; 
PLANES{4} = [4 5] ;

PLANES{5} = [1 6] ;
PLANES{6} = [2 6] ;
PLANES{7} = [3 6] ; 
PLANES{8} = [4 6] ;

PLANES{9} =  [1 2 ] ;
PLANES{10} = [1 4 ] ;
PLANES{11} = [2 3 ] ;
PLANES{12} = [3 4 ] ;



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