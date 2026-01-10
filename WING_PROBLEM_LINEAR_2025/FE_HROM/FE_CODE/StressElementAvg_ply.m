function  stressPLY = StressElementAvg_ply(DATAOUT)
% COmpute average element stress for the elements of each ply 
nstress = 6 ; 
ngaus =  length(DATAOUT.stress)/nstress;  
nelem = length(DATAOUT.MaterialType) ; 
ngausE = ngaus/nelem ; 
%%% AVERAGE STRESSES ON ELEMENTS (GIVEN STRESSES AT ALL GAUSS POINTS)
stressEL = Stress_Elements_Average(DATAOUT.wSTs,DATAOUT.stress,ngausE,nstress,nelem) ; 
%%%%
% Average separated for different plies 
nmat = length(unique(DATAOUT.MaterialType) ) ; 

stressPLY ={} ; 
for imat = 1:nmat 
    elem = find(DATAOUT.MaterialType==imat) ; 
    stressPLY{imat} = stressEL(:,elem) ;
end
 