
function AuxDATA=WriteAuxFdNamesNEW(FdnamesInputs,varginOR,DATAPath,varargin)
% It writes an script for reading inputs
% Cruciall --> For reading inputs
% Example
% WriteAuxFdNames(DATAPath,'AuxFile',FdnamesInputs) 
% AuxFile ; --> Scrip
% see PlotVarProj
% 
% if nargin == 3
%     DATAPath = getenv('matpathdata');
% end

% Default_Inputs  
NameStrInputs = 'StrInputs';
NameOk        = 'okInputs' ;
DATAPath      = getenv('matpathdata');
PathGidCode = getenv('GID_CODE_M');
    
% End Default Inputs
vararginK = varargin;
[OKNameStrInputs,NameStrInputs,vararginK]=FindInArgOUT(vararginK,'NameStrInputs',NameStrInputs);
[OKNameOk,NameOk,vararginK]=FindInArgOUT(vararginK,'NameOk',NameOk);
[OKk,DATAPath,vararginK]=FindInArgOUT(vararginK,'DATAPath',DATAPath);

%fid = fopen(AuxFdNames_file,'w');
str = '';
for i = 1:length(FdnamesInputs) 
    str = [str,FdnamesInputs{i},','];
end
AuxDATA{1} =['DefInputs = {',str(1:end-1),'} ;'];
%eval(str) ;
%ppp = fprintf(fid,'%s\n',str);
% Default inputs
% we change those variables appearing in varargin
str2 = ['[',NameStrInputs,' ',NameOk,']','  = VarArginInput(DefInputs,FdnamesInputs,varginOR);'];
AuxDATA{2}=str2;
%ppp = fprintf(fid,'%s\n',str2);
% Now we write something like:
%
% StrInputs.Var = Var ...
%
%
for i = 1:length(FdnamesInputs)     
    str = [FdnamesInputs{i},'=','StrInputs.',FdnamesInputs{i},' ;'];
    AuxDATA{2+i}=str;
    %ppp = fprintf(fid,'%s\n',str);
end

%fod = fclose(fid);