function CheckMeshJoint(NAMEMESH_SLICE,itype,CNint,MaterialType1D,DATAIN,SLICEglo,rotMATglo,coorINT)
% Check whether NAMEMESH_SLICE{itype}  exists. If not, it helps the user to
% create the geometry
% JAHO, 7-Jan-2018

if exist(NAMEMESH_SLICE{itype},'file') ==0
    warning(['Joint ', NAMEMESH_SLICE{itype} ,' has not been yet created'])
    ok = menu(['Joint ', NAMEMESH_SLICE{itype} ,' does not exist. Do you want to create it ? '],['YES'],['NO'])
    %% Creating new joint
    if ok==1
        CreateNewJoint(NAMEMESH_SLICE,CNint,MaterialType1D,itype,DATAIN,SLICEglo,rotMATglo,coorINT)
    else
        if   exist([NAMEMESH_SLICE{itype}(1:end-4),'.gid'],'file') >0
            TTT = ['gid ',NAMEMESH_SLICE{itype}(1:end-4),'.gid'] ;
            unix(TTT) ;
            error('Edit the GID file, export the mesh and run again')
        else
            error('Non-existing joint')
        end
    end
    
    
    
end