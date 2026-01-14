function  GidMesh2DFE(NameFile,COOR,CONNECT,NAMEPROJ,MaterialType,TypeElement)

%%% INPUTS %%%%%%
numer_nodesTOT = 1:size(COOR,1) ; 
% TypeElement='Quadrilateral' ; 
% %% EXTRACTING INPUTS
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varginOR = varargin ;
% FdnamesInputs = {'numer_nodes','TypeElement'};
% AuxDATA=WriteAuxFdNamesNEW(FdnamesInputs,varginOR);
% for id = 1:length(AuxDATA);
%     eval(AuxDATA{id});
% end
%E_JAHO

% try
% DATA = DefaultField(DATA,'NODES',1:size(COOR,1) ) ; 
% catch
%     DATA.NODES = 1:nnode  ; 
% end
%numer_nodes =   DATA.NODES  ; 

elem_type = TypeElement;
 NNode=size(CONNECT,2);
 ndime = size(COOR,2) ;
 nnod = size(COOR,1) ;
 npe =  size(CONNECT,2) ;
 nElem = size(CONNECT,1) ;
 
 [ppp,~,~] = fileparts(NameFile) ; 
 if ~exist(ppp,'file')
     mkdir(ppp); 
 end
 
 
fid = fopen(NameFile,'wt');

% 1.- Header with 6 free lines
[dummy, NAMEPROJ ]= fileparts(NAMEPROJ) ; 
fprintf(fid,' # ================================================== \n',[]);
fprintf(fid,[' # \n'],[]);
fprintf(fid,' # POSTPROCESSING WITH GID - MESH FILE \n',[]);
fprintf(fid,[' # EXAMPLE NAME: ' NAMEPROJ ' \n'],[]);
fprintf(fid,' # ================================================== \n',[]);
fprintf(fid,[' MESH "' NAMEPROJ '" dimension ' num2str(ndime) ' Elemtype ' elem_type ...
              ' Nnode ' num2str(NNode) '\n'],[]);
fprintf(fid,' Coordinates \n',[]);
fprintf(fid,' # node number   coordinate_x  coordinate_y  coordinate_z \n',[]);

% 4.- Coordinates 
format_xnod = [];
for k=1:ndime,
   format_xnod = [format_xnod, '%15.5e'];
end
format_xnod=[' %10i ' format_xnod, '\n'];
fprintf(fid,format_xnod,[numer_nodesTOT',COOR(:,1:ndime) ]');

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

fprintf(fid,format_icone,[(1:nElem)',CONNECT, MaterialType]');

fprintf(fid,' end elements ',[]);

fclose(fid);