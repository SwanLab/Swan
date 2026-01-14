function [OkFindRes ,nrow] = ReadUntCount(fid,Word,ncolumn)
% Given a file "fid", it reads until find the word "Word" 
%  For example 
%  6 7 8
%  5 6 7
%  3 4 5
%  end
%  Word = end , ncolumn = 3 --> nrow = 3
%

OkFindRes = 0;
nrow = 0;

while OkFindRes==0
    leido=fscanf(fid,'%s',1);
    if strcmpi(Word,leido);
        % Fin del bucle
        OkFindRes=1;
    else
        dummy = fscanf(fid,'%s',ncolumn-1) ;
        nrow=nrow+1;
      
        
    end
end
