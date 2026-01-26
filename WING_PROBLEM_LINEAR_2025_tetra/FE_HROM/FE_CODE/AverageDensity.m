function   DATAOUT =   AverageDensity(DATAOUT,MATERIAL,DATA) ; 

 % ---------------------------------
 % Computing  AVERAGE DENSITY  
 % ---------------------------------   
%dbstop('7')
 if nargin ==0
     load('tmp.mat')
 end
%  
nstress = 6 ; 
ngaus =  length(DATAOUT.stress)/nstress;  
nelem = length(DATAOUT.MaterialType) ; 
ngausE = ngaus/nelem ; 
% volume of each element 
weig = DATAOUT.wSTs ;
volE = reshape(weig,ngausE,[]) ; 
volE = sum(volE) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
elemFIB = find(DATAOUT.MaterialType == MATERIAL.FIBER.INDEX) ;  % Fiber elements
elemMAT= find(DATAOUT.MaterialType == MATERIAL.MATRIX.INDEX) ;  % Matrix elements 


volFIB = sum(volE(elemFIB)) ; 
volMAT =  sum(volE(elemMAT)) ; 

volTOT = volFIB+volMAT ; 

fracVOL.FIBER = volFIB/volTOT ; 
fracVOL.MATRIX = volMAT/volTOT ; 

densCOMP = (volFIB*MATERIAL.FIBER.DENSITY + volMAT*MATERIAL.MATRIX.DENSITY)/volTOT ; 

DATAOUT.densCOMP = densCOMP ;
