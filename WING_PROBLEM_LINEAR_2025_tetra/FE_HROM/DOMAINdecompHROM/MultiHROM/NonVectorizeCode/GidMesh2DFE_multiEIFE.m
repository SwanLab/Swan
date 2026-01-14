function  IND_ELEM_MESHES = GidMesh2DFE_multiEIFE(NameFile,COOR,CONNECT_glo,NAMEPROJ,MaterialType,TypeElement,NAMEMESH)
if nargin == 0
    load('tmp1.mat')
end

%%% INPUTS %%%%%%
numer_nodesTOT = 1:size(COOR,1) ;


fid = fopen(NameFile,'wt');
fprintf(fid,' # ================================================== \n',[]);
fprintf(fid,[' # \n'],[]);
fprintf(fid,' # POSTPROCESSING WITH GID - MESH FILE \n',[]);
fprintf(fid,[' # EXAMPLE NAME: ' NAMEPROJ ' \n'],[]);
fprintf(fid,' # ================================================== \n',[]);
iacumELEM = 0 ; 
IND_ELEM_MESHES ={} ;
for imesh = 1:length(CONNECT_glo)
    
    elem_type = TypeElement{imesh};
    CONNECT= CONNECT_glo{imesh};
    NNode=size(CONNECT,2);
    ndime = size(COOR,2) ;
    nnod = size(COOR,1) ;
    npe =  size(CONNECT,2) ;
    nElem = size(CONNECT,1) ;
    NAMEMESH_loc = NAMEMESH{imesh} ;
    MaterialType_loc = MaterialType{imesh} ;
    IND_ELEM_MESHES{imesh} = [iacumELEM+1:iacumELEM+nElem] ; 
    % 1.- Header with 6 free lines
    
    
    fprintf(fid,[' MESH "' NAMEMESH_loc '" dimension ' num2str(ndime) ' Elemtype ' elem_type ...
        ' Nnode ' num2str(NNode) '\n'],[]);
    fprintf(fid,' Coordinates \n',[]);
    if imesh ==1
        fprintf(fid,' # node number   coordinate_x  coordinate_y  coordinate_z \n',[]);
        
        % 4.- Coordinates
        format_xnod = [];
        for k=1:ndime,
            format_xnod = [format_xnod, '%15.5e'];
        end
        format_xnod=[' %10i ' format_xnod, '\n'];
        fprintf(fid,format_xnod,[numer_nodesTOT',COOR(:,1:ndime) ]');
    end
    fprintf(fid,' end coordinates \n',[]);
    
    fprintf(fid,' Elements \n',[]);
    
    format_title = ' # element ';
    for k=1:npe,
        format_title=[format_title, ' node_' num2str(k)];
    end
    format_title=[format_title, ' material number \n'];
    
    fprintf(fid,format_title,[]);
    
    format_icone = [];
    for k=1:npe,
        format_icone=[format_icone, '%10i'];
    end
    format_icone=['%10i ' format_icone, '     %10i ' '\n'];
    
    fprintf(fid,format_icone,[((1:nElem)+iacumELEM)',CONNECT, MaterialType_loc]');
    
    fprintf(fid,' end elements \n',[]);
    iacumELEM = iacumELEM + nElem ; 
end

fclose(fid);