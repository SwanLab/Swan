function DATA = ClusterComputeSize(LIMIT,DATA)


%nclusters = 2; 
% Assume we have 5 Gb available. Matrices over, say, 100 Mb should be
% written to disk. This is the limit
%LIMIT = DATA.LIMIT_mbytes_matrices ; % ;Mb
% What is the size of nodal variables
%
nsizeNOD = DATA.MESH.ndof*length(DATA.STEPS)*8*1e-6 ;
nsizeGAUSS = DATA.MESH.ndofSTRESS*length(DATA.STEPS)*8*1e-6 ;

nclusters = max(ceil(nsizeNOD/LIMIT),ceil(nsizeGAUSS/LIMIT)) ;
DATA.nclusters = nclusters; 

DATA.STORE.NSTEPS_CLUSTER = cell(1,nclusters) ;
% We have to divide  1:length(DATA.STEPS) into nclusters
FREQ  = ceil(length(DATA.STEPS)/nclusters) ;
iini = 1;
%ifin = FREQ ;
for i=1:nclusters
    ifin = iini + FREQ-1 ;
    ifin = min(ifin,length(DATA.STEPS)) ;
    DATA.STORE.NSTEPS_CLUSTER{i} = iini:ifin ;
    iini = ifin +1 ;
end