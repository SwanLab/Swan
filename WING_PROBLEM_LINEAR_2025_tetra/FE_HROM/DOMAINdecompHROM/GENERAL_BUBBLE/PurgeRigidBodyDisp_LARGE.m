function SNAPcompl = PurgeRigidBodyDisp_LARGE(DATAoffline,SNAPcompl,Mdom,MESH)
% Purging rigid body displacements from displacement snapshots
% Large rotations case
% Theory based on
% /home/joaquin/Desktop/CURRENT_TASKS/PAPERS_2020_onwards/12_EIFEM_EXTENS/EIFEM_largeROTfinal.pdf
% JAH0, 18-May-2025, Sunday, Secrets by Farga, Barcelona.
% ------------------------------------------------------------------------------------

if nargin == 0
    load('tmp.mat')
end
 


 
for iproj = 1:length(SNAPcompl)
    
    % PURGE ROTATIONS
    for isnap = 1:size(SNAPcompl{iproj},2)
        U = SNAPcompl{iproj}(:,isnap);
        U = reshape(U,size(MESH.COOR,2),[])' ;  
        [Q, t, X_rigid, Udeform] = ExtractRigidMotion_Vectorized(Mdom,  MESH.COOR, U) ;
        Udeform = Udeform' ; 
        SNAPcompl{iproj}(:,isnap) = Udeform(:) ; 
    end
    
    
end
