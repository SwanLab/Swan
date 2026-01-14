function [  TypeOfSlice] = WhichTypeOfSlices(MESH1D)

% Type of Slice  within each beam
% -----------------
% ------------------------
% Loop over entities loop
% --------------------------
if nargin == 0
    load('tmp1.mat')
end

[nelem nnodeE]= size(MESH1D.CN) ;
TypeOfSlice = zeros(nelem,1) ;  % zero by default
for iENTITY  = 1: length(MESH1D.PROP)
    %   switch MESH1D.PROP(iENTITY).TYPE ;
    %      case {'BEAM'}
    itypeGLO = MESH1D.PROP(iENTITY).INDEX_ENTITY ;
    % Default properties for this type of beam
    PROP = MESH1D.PROP(iENTITY) ;
    % Elements pertaining to this type of beam
    ELEMS = find(MESH1D.MaterialType == iENTITY) ;
    
    if length(itypeGLO) == 1
        TypeOfSlice(ELEMS)  = itypeGLO  ;
        
    else
        % MIXED BEAMS
        % ----------
        %%% SPECIFIC ELEMENTS of type slice "itypeSLICE"
        %MESH1D.PROP(itype).ELEMENT_OF_EACH_ENTITY = {1:3:114,2:3:114,3:3:114};
        PROP = DefaultField(PROP,'ELEMENT_OF_EACH_ENTITY',[]) ;
        for itypeSLICE = 1:length(itypeGLO)
            if isempty(PROP.ELEMENT_OF_EACH_ENTITY) % || isempty(PROP.ELEMENT_OF_EACH_ENTITY{itypeSLICE})
                error('Specify the 1D elements associated with each slice !!! ')
            else
                if isempty(PROP.ELEMENT_OF_EACH_ENTITY{itypeSLICE})
                    ELEM_LOC = ELEMS ;
                else
                    ELEM_LOC = PROP.ELEMENT_OF_EACH_ENTITY{itypeSLICE} ;
                end
                
                TypeOfSlice(ELEM_LOC) = itypeGLO(itypeSLICE) ;
            end
        end
    end
    
    %  end
end

end
