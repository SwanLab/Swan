clc
clear all

nvol = 10 ;

TXT ={} ;

for ivol = 1:nvol
    
    TXT{end+1}  = ['*Set Cond LINE-',num2str(ivol),' *elems'] ;
    TXT{end+1}  = '*loop elems *OnlyInCond' ;
    TXT{end+1}  =  '*format "%5i%5i"' ;
    TXT{end+1} ='*elemsNum *cond(1,int)' ;
    TXT{end+1} ='*end' ;
    TXT{end+1} ='   ' ;
end

NAME = [cd,filesep,'tmp.txt'];

fid =fopen(NAME,'w');
for i = 1:length(TXT)
    
    fprintf(fid,'%s\n',[TXT{i}]);
    
end
fod =fclose(fid);

open(NAME)