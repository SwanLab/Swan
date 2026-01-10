function [  Elements2Print] = WhichBeamElementsToPrint(MESH1D)

% Slice domains to be printed 
% -----------------
% ------------------------
% Loop over entities loop
% --------------------------
if nargin == 0
    load('tmp1.mat')
end

[nnode ndim]= size(MESH1D.COOR) ;
[nelem nnodeE]= size(MESH1D.CN) ;
Elements2Print = [] ;
 

for iENTITY  = 1: length(MESH1D.PROP)
    itypeGLO = MESH1D.PROP(iENTITY).INDEX_ENTITY ;
    switch MESH1D.PROP(iENTITY).TYPE ;
        case {'BEAM'}
            % Default properties for this type of beam
            PROP = MESH1D.PROP(iENTITY) ;
            IND_BEAMS = find(MESH1D.MaterialType == iENTITY) ;
            PROP = DefaultField(PROP ,'DOMAINSprint',[]) ;
            % Selection method: default: Specify % of slices to be printed (10%, for instance)
            PROP.DOMAINSprint = DefaultField(PROP.DOMAINSprint ,'TYPE_SELECTION','INDICES_SLICES') ;
            % Other option: INDICES_SLICES --> This should be specify by the
            % user. Global indices read from the 1D mesh
            PROP.DOMAINSprint = DefaultField(PROP.DOMAINSprint ,'VARIABLES',IND_BEAMS) ; % Variables defininig selection
            
            switch   PROP.DOMAINSprint.TYPE_SELECTION
                case 'INDICES_SLICES'
                    Elements2PrintLocal =  PROP.DOMAINSprint.VARIABLES ;
                case 'PERCENTAGE'
                    % List of separate portions for this type of beam
                    ELEMS_each_portion = MESH1D.INFOLINES.ELEMENTS{iENTITY}   ;
                    Elements2PrintLocal = [] ;
                    % Loop over portion of beams
                    for iportion =  1:length(ELEMS_each_portion)  % Loop over portions of beams
                        ELEMS = ELEMS_each_portion{iportion} ; % List of 1D elements corresponding to a given portion
                        Percentage=PROP.DOMAINSprint.VARIABLES;
                        
                        if isempty(Percentage) | Percentage >=100
                            ToPrint = ELEMS ;
                        else
                            nprint = floor(length(ELEMS)*Percentage/100) ;
                            FREQ = ceil(length(ELEMS)/nprint) ; 
                            PRINT_ELEMENTS = 1:FREQ:length(ELEMS) ;
                            ToPrint = ELEMS(PRINT_ELEMENTS) ;
                            % Final and initial domains are always printed
                            % (to check that connection with joints are correct)
                            ToPrint = unique([ToPrint;ELEMS(1);ELEMS(end)]) ;
                            
                        end
                        
                        Elements2PrintLocal = [Elements2PrintLocal;ToPrint ] ;
                    end
                    
                otherwise
                    error('Option not implemented')
            end
            
            Elements2Print = [Elements2Print ;Elements2PrintLocal ] ;
            %   TypeElements2Print = [TypeElements2Print ; iENTITY*ones(size(Elements2PrintLocal))] ;
    end
end

%Elements2Print = [Elements2Print,TypeElements2Print]  ;


%
%    for jslice = 1:length(itypeGLO)
%                 itypeSLICE = itypeGLO(jslice) ;
%
%
%
%                 % We select here the elements whose normals will be printed
%
%
%
%
%
%
%                 Elements2Print  = [Elements2Print ; Elements2PrintLoc] ;
%                 SliceElementToPrint = [SliceElementToPrint ; itypeSLICE*ones(size(Elements2PrintLoc))] ;
%                 TypeElements2Print = [TypeElements2Print ; iENTITY*ones(size(Elements2PrintLoc))] ;
%    end
%
% Elements2Print = [Elements2Print,TypeElements2Print , SliceElementToPrint] ;