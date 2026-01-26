function CheckMeshJoint(itypeJOINT,MESH1D,MESH3D,DATAIN)
% Check whether NAMEMESH_SLICE{itype}  exists. If not, it helps the user to
% create the geomtry
% J AHO, 10-Jan-2018
if nargin == 0
    load('tmp1.mat')
end

NAMEmsh = [MESH3D.JOINTS(itypeJOINT).NAME,'.msh'] ; 
if exist(NAMEmsh,'file') ==0
    warning(['Joint ', MESH3D.JOINTS(itypeJOINT).NAME,' has not been yet created'])
    ok = menu(['Joint ', MESH3D.JOINTS(itypeJOINT).NAME,' does not exist. Do you want to create it ? '],['YES'],['NO'])
    %% Creating new joint
    if ok==1
        WRITE_BATCH_FILE = 1; 
        CreateNewJoint(MESH1D,MESH3D,itypeJOINT,WRITE_BATCH_FILE) ;
    else
        if   exist([JOINT.NAME(1:end-4),'.gid'],'file') >0
            TTT = ['gid ',JOINT.NAME(1:end-4),'.gid'] ;
            unix(TTT) ;
            error('Edit the GID file, export the mesh and calculation file and run again')
        else
            error('Non-existing joint')
        end
    end
    
    
    
end