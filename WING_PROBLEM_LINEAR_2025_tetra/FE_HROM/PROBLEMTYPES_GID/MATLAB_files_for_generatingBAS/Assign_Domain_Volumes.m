clc
clear all

nvol = 48 ;

TXT ={} ;

 
 

for ivol = 1:nvol
    
    TXT{end+1}  = ['Mescape'] ;
    TXT{end+1}  = 'escape escape escape escape' ;
    TXT{end+1}  =  ['Mescape Data Conditions AssignCond DOMAIN-',num2str(ivol),' Change ',num2str(ivol)] ;
    TXT{end+1} =num2str(ivol) ;
    TXT{end+1}  =   'escape';
  
 end

NAME = [cd,filesep,'tmp.txt'];

fid =fopen(NAME,'w');
for i = 1:length(TXT)
    
    fprintf(fid,'%s\n',[TXT{i}]);
    
end
fod =fclose(fid);

open(NAME)