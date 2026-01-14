function [OkFindRes ,RemLine] = ReadUntilToken(fid,Token,N)
% Given a file "fid", it reads until find the word "TOKEN" (or until read N lines)
if nargin == 2
    N = 1E20;
end


TokenRead = 'dummy';
nread = 0;
while strcmpi(TokenRead,Token)~=1 & (TokenRead~=-1 ) & nread < N
    TokenRead = fgets(fid);
    nread = nread +1;
    [TokenRead RemLine] = strtok(TokenRead);
    if isempty(TokenRead)
        TokenRead = 'pasa_al_siguienteee';
    end
end

if strcmpi(TokenRead,Token)== 1
    OkFindRes = 1;
else
    OkFindRes = 0;
end
