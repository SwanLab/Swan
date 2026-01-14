function [NshapeFINAL,COORnodes,OTHER_OUTPUT] =  IntfModes_AUXalex_EIFEM(INFO_INTERFACE_MESH,CENTROID,COORbnd)
% JAHO, 8-MAY-2025, THURSDAY, upc TERRASSA
% Computation of Nshape for the interface ficti. modes, in EIFEM
% mIXED arbitrary interolation/LINEAR INTERPOLATIONS, 2D.
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/109_EIFEM_largeROT/12_REPETITIVE_TRAIN.mlx
% --------------------------------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end
OTHER_OUTPUT = [] ;

INFO_INTERFACE_MESH = DefaultField(INFO_INTERFACE_MESH,'xi_nodes',[]) ;
xi_nodes = INFO_INTERFACE_MESH.xi_nodes;



COORnodes = INFO_INTERFACE_MESH.COOR; %(:,2:end) ;
ElemBnd= INFO_INTERFACE_MESH.LINES ;
nNODES_with_repetition = sum(cellfun(@numel, ElemBnd));

nmodes = size(COORnodes,1) ; % Number of modes
npoints = size(COORbnd,1) ;



Nshape = zeros(npoints,nNODES_with_repetition) ;

% cOORDINATES REFERRED TO THE CENTROID
ndim = size(CENTROID,2);
for idim = 1:ndim
    COORnodes(:,idim) = COORnodes(:,idim) - CENTROID(idim) ;
end

xiEDGES = cell(size(ElemBnd)) ;
IndPointsBNDedge = cell(size(ElemBnd)) ;
NumberNodesPerEdge = zeros(size(ElemBnd)) ;


for ielem = 1:length(ElemBnd)
    
    
    
    xi_nodes = linspace(-1,+1,length(ElemBnd{ielem})) ;
    
    [Nshape,xiEDGES{ielem},IndPointsBNDedge{ielem}] =   QuarticShapeFun5nodes(ElemBnd{ielem},COORnodes,COORbnd,Nshape,xi_nodes) ;
    
    NumberNodesPerEdge(ielem) = length(ElemBnd{ielem});
    
    
    
end
OTHER_OUTPUT.xiEDGES = xiEDGES ;
OTHER_OUTPUT.IndPointsBNDedge = IndPointsBNDedge ;
OTHER_OUTPUT.NumberNodesPerEdge = NumberNodesPerEdge ;

% Now we have to transform Nshape into NshapeFINAL
% How to do it ? I'd construct with a transformation marix
TransformMatr = zeros(nNODES_with_repetition,nmodes) ;
% This matrix is equal to the identiy in case of disjoint line segments
% By contrast, if lines are connected, the contribution of the share node
% should be weighted by 0.5 (assuming that all nodes are numbered consecutively)
% So we have to transvers all nNODES_with_repetition and choose the
% corresponding weights
iNEW = 1;
iOLD = 1 ;

if ElemBnd{end}(end) == ElemBnd{1}(1)
    IS_repeated = 1;
    TransformMatr(end,1) = 0.5 ;
end

for  iLINES = 1:length(ElemBnd)
    % Loop over number of segments
    LineLoc = ElemBnd{iLINES} ;
    for inodeLOC = 1:length(LineLoc)
        % Loop over nodes of a given segment
        inodeGLO_current = LineLoc(inodeLOC) ;
        % First we have to check whether this
        % node is shared with another segment
        % This can only occur if it is either the
        % first node or last node of the segment.
        if inodeLOC ==1
            if IS_repeated == 1
                % This has been checked when inspecting the previous
                % segement in the sequence
                TransformMatr(iNEW,iOLD) = 0.5 ;
            else
                TransformMatr(iNEW,iOLD) =1;
            end
            iOLD = iOLD + 1;
        elseif inodeLOC == length(LineLoc)
            % Check which is the next node in the sequence
            if iLINES < length(ElemBnd)
                inodeGLO_next = ElemBnd{iLINES+1}(1) ;
                if inodeGLO_next == inodeGLO_current
                    IS_repeated = 1 ;
                    TransformMatr(iNEW,iOLD) = 0.5;
                else
                    IS_repeated = 0 ;
                    TransformMatr(iNEW,iOLD) = 1;
                    iOLD = iOLD + 1;
                end
            else
                % This is the last node, it has been already assigned
            end
        else
            TransformMatr(iNEW,iOLD) = 1;
            iOLD = iOLD + 1;
        end
         iNEW = iNEW +1;
    end
   
end


NshapeFINAL = Nshape*TransformMatr ;
