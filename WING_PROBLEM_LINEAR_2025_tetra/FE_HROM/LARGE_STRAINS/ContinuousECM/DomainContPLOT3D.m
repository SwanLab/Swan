function DomainContPLOT3D(AUXVAR,MESH,VAR_SMOOTH_FE)

if nargin == 0
    load('tmp2.mat')
   figure(1)
     hold on
     axis equal
%     AUXVAR.DATALOC.SHOW_MESH_DOMAIN  = 1; 
 
end


AUXVAR.DATALOC = DefaultField(AUXVAR.DATALOC,'ContourMeshDomain',[]) ;
meshDOM = AUXVAR.DATALOC.ContourMeshDomain ;
AUXVAR.DATALOC = DefaultField(AUXVAR.DATALOC,'SHOW_MESH_DOMAIN',0) ;

 
if ~isempty(meshDOM) && AUXVAR.DATALOC.SHOW_MESH_DOMAIN == 0
    LinesToPlot = [
        11     1     9    22    23    10    21    27    11     1     9    22
        1     9    22    11    10    21    27    23    23    10    21    27]' ; % This is for quadratic elements
    meshDOM = DefaultField(meshDOM,'LinesToPlot',LinesToPlot) ;
     
    if size(meshDOM.LinesToPlot,1) == 2
        meshDOM.LinesToPlot =  meshDOM.LinesToPlot' ; 
    end
    
    for iplot = 1:size(meshDOM.LinesToPlot,1)
        XX = cell(1,3) ;
        for idim = 1:3
            XX{idim} = [meshDOM.COOR( meshDOM.LinesToPlot(iplot,1),idim),meshDOM.COOR(meshDOM.LinesToPlot(iplot,2),idim)] ;
        end
        plot3(XX{1},XX{2},XX{3},'LineWidth',0.5,'Color',[0 0 0] )
%          if iplot == 1
%              text(XX{1}(1),XX{2}(1),XX{3}(1),num2str(meshDOM.LinesToPlot(iplot,1)))
%              text(XX{1}(2),XX{2}(2),XX{3}(2),num2str(meshDOM.LinesToPlot(iplot,2)))
%          else
%               text(XX{1}(2),XX{2}(2),XX{3}(2),num2str(meshDOM.LinesToPlot(iplot,2)))
%          end
    end
    
elseif AUXVAR.DATALOC.SHOW_MESH_DOMAIN == 1 
    error('This option does not work properly')
    
for elemLOC = 1:size(VAR_SMOOTH_FE.CN,1)  
 %   elemLOC = INDused(ielem) ;
    NODESloc  = VAR_SMOOTH_FE.CN(elemLOC,VAR_SMOOTH_FE.IND_POLYG_ELEMENT)' ;
    COORloc = VAR_SMOOTH_FE.COOR(NODESloc,:) ;
    plot3(COORloc(:,1),COORloc(:,2),COORloc(:,3),'Color',[0.1,0.1,0.1]) ;
end
    
     
end
