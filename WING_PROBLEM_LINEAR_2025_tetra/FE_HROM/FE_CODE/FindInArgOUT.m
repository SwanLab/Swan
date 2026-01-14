function [OKFIELD ,VALUE,vararginOUT]=FindInArgOUT(vararginCELL,Field,DefaultValue)
% Similar to FindInArg, bu it optimizes the process.
%   [OKFIELD VALUE]=FindInArg(vararginCELL,Field)
%   For using "varargin". It looks for string "Field" in vararginCELL and
%   returns if it has been found and it's value (if it exists)
%

if nargin == 0
    vararginCELL ={[0 0],'Ndiv',20,'Ellipse','Radius',[2 3]};
    Field = 'Ndiv';
    addpathGLO;
elseif nargin == 2
    DefaultValue = 0;
end

[OKFIELD,nsortFIELD]=FndStrInCell(vararginCELL,Field);

if OKFIELD == 1
    VALUE = vararginCELL{nsortFIELD+1};
    % Reshaping 
    vararginCELL(nsortFIELD:nsortFIELD+1) = [];  
else
    VALUE = DefaultValue;   
end
 vararginOUT = vararginCELL;
