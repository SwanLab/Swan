function [StrInputs, okInputs] = VarArginInput(DefInputs,FdnamesInputs,vararginOR)
% It finds variables in vararginOR and builds two structures StrInputs,okInputs
% Example 
% vararginOR = {'a',5};
% DefInputs = {'a','b','c'}
% FdnamesInputs = {'letra1','letra2','letra3'}
%StrInputs =

% letra1: 'alpha'
% letra2: 'b'
% letra3: 'c'
% 
% 
% okInputs =
% 
% letra1: 1
% letra2: 0
% letra3: 0

%

if nargin == 0
    vararginOR = {'letra1','alpha'};
    DefInputs = {'a','b','c'};
     FdnamesInputs = {'letra1','letra2','letra3'};
%
end
StrInputs = cell2struct(DefInputs,FdnamesInputs,2);
okInputs  = StrInputs; 

% Now we replace the non-defaults inputs 
for ifield = 1:length(FdnamesInputs)    
    [ok,NewInp,vararginOR]=FindInArgOUT(vararginOR,FdnamesInputs{ifield},DefInputs{ifield});
    StrInputs=setfield(StrInputs,FdnamesInputs{ifield},NewInp);
    okInputs=setfield(okInputs,FdnamesInputs{ifield},ok);
end
