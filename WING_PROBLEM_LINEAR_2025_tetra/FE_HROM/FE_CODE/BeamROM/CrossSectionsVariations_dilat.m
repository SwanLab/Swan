function   [A,D,a0] = CrossSectionsVariations_dilat(MESH1D,MESH3D,DATAIN)

% TRANSFORMATION MATRICES, See BeamROM.pdf
% JAHO, 23-January-2018
% -----------------
% ------------------------
% Loop over entities loop
% --------------------------
if nargin == 0
    load('tmp1.mat')
end


[nnode ndim]= size(MESH1D.COOR) ;
[nelem nnodeE]= size(MESH1D.CN) ;
% Parameters defining transformation -> x = a0 + A*X + X1*D*X, see
% BeamROM.pdf
% ----------------------------------------------------------------
a0 = zeros(ndim,nelem) ;
A  = zeros(ndim,nelem*ndim) ;  % Dilatational matrices
D = zeros(ndim,nelem*ndim)  ;
for imat1D  = 1: length(MESH1D.PROP)
    itypeGLO = MESH1D.PROP(imat1D).INDEX_ENTITY ;
    if ~isempty(itypeGLO)
        switch MESH1D.PROP(imat1D).TYPE ;
            case 'BEAM'
                % Loop over "subentities" (beams)
                % -----------------------
                nbeams = length(MESH1D.INFOLINES.ELEMENTS{imat1D}) ;  % Number of "beams" of this type
                
                
                for ibeam = 1:nbeams
                    [A,D,a0] = TransformationMatrices(MESH1D,MESH3D,imat1D,itypeGLO,A,D,a0,ibeam) ;
                end
        end
    end
end
