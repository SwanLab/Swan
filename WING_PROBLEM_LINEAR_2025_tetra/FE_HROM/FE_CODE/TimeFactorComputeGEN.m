function [TIME_FACTOR,STEPS_TO_STORE,STEPS_TO_PRINT,TIME_DISCRETIZATION] ...
    = TimeFactorComputeGEN(FACTORS,NSTEPS,NSTEPS_STORE_input,NSTEPS_PRINT_input)

TIME_FACTOR = {} ; 
TIME_ACUM =  0  ; 
NSTEPS_STORE_SUM = 0 ; 
NSTEPS_STORE = {}; 
NSTEPS_PRINT = {} ; 

 

for itime = 1:length(NSTEPS)
    if itime ==1
        MAS = 0 ; 
    else
        MAS  =1; 
    end
    TIME_FACTOR{itime} =   linspace(FACTORS(itime),FACTORS(itime+1),NSTEPS(itime)+MAS) ; 
    NSTEPS_STORE{itime} = NSTEPS_STORE_SUM + ceil(linspace(1,NSTEPS(itime),NSTEPS_STORE_input(itime))) ; 
    NSTEPS_STORE_SUM = NSTEPS_STORE_SUM + NSTEPS(itime) ; 
end
 
TIME_FACTOR = (cell2mat(TIME_FACTOR)) ; 
BBB = diff(TIME_FACTOR) ; 
XXX = find(BBB==0) ; 
TIME_FACTOR(XXX+1) = []  ; 

TIME_DISCRETIZATION = abs(diff(TIME_FACTOR)) ; 
TIME_DISCRETIZATION = [0 cumsum(TIME_DISCRETIZATION)] ; 

 
 
[~,NSTEPS_STORE_input] = cellfun(@size,NSTEPS_STORE) ; 
PRINT_A = 0 ; 
 for itime = 1:length(NSTEPS_STORE_input)
    NSTEPS_PRINTloc =  ceil(linspace(1,length(NSTEPS_STORE{itime}),NSTEPS_PRINT_input(itime))) ; 
    NSTEPS_PRINT{itime} =NSTEPS_PRINTloc + PRINT_A ;
    PRINT_A = PRINT_A +    NSTEPS_STORE_input(itime)  ; 
 end

 STEPS_TO_PRINT=  unique(cell2mat(NSTEPS_PRINT)) ; 
 STEPS_TO_STORE=  unique(cell2mat(NSTEPS_STORE)) ; 