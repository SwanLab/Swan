clc
clear all

nvol = 50 ;
nini = 50 ; 

TXT ={} ;

for ivol = 1:nvol
    
    TXT{end+1}  = ['NUMBER: ',num2str(nini+ivol),' CONDITION: LINE-',num2str(ivol)] ;
    TXT{end+1}  = 'CONDTYPE: over lines' ;
  %  TXT{end+1}  =  'CONDMESHTYPE: over face elements' ;
   TXT{end+1}  =  'CONDMESHTYPE: over nodes' ;
    TXT{end+1} ='CANREPEAT: yes' ;
    TXT{end+1} ='QUESTION: ASSOCIATED_INDEX' ;
    TXT{end+1} =['VALUE: ',num2str(ivol)] ;
    TXT{end+1} ='END CONDITION' ;
    
end

NAME = [cd,filesep,'tmp.txt'];

fid =fopen(NAME,'w');
for i = 1:length(TXT)
    
    fprintf(fid,'%s\n',[TXT{i}]);
    
end
fod =fclose(fid);

open(NAME)