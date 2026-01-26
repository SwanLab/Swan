function [OK,nsort]=FndStrInCell(InpCell,InpStr)
%        [OK,nsort]=FndStrInCell(InpCell,InpStr)
% It finds string "InpStr" in the string cell InpCell ;
%profile on 

if nargin == 0
    InpCell = {'ABC','DEF','GDG'};
    InpStr  = 'DEF';
end 
    

vector = strcmp(InpCell,InpStr) ;
nsort = find(vector == 1); 
if isempty(nsort)
    OK = 0 ; nsort = 1 ; 
else
    OK = 1 ;
    nsort =  nsort(1) ;
end



% 
% InpCell_i = 'kdkdkddkdkdk';
% ncell = length(InpCell);
% istr = 1;
% nsort = 1;
% while strcmp(InpCell_i,InpStr)==0 & istr<=ncell
%     InpCell_i = InpCell(istr);
%     istr = istr + 1;
% end
% 
% 
% OK = strcmp(InpCell_i,InpStr) ;
% 
% if OK == 0
%     nsort = 1;
% else
%     nsort = istr-1;
% end

%profile report