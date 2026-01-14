function  [VARNODES,VARGAUSS ]= GaussToNodesSmooth(VARGAUSS,Nst,DATAINloc) ;
%  TRANSPORTING  VARIABLES FROM GAUSS POINTS TO NODES
% JAHO, 17-Apr-2020  (35th of confinement because of COVID19)
% ----------------------------------------------------------------
if nargin == 0
    load('tmp1.mat')
end

DATAINloc = DefaultField(DATAINloc,'GaussPointsToNodesMethod','GlobalSmoothinMassMatrix') ;

ncomponent = DATAINloc.ncomponent ;
nnodes = size(Nst,2) ;
VARNODES = zeros(nnodes*ncomponent,size(VARGAUSS,2)) ;


DATAINloc = DefaultField(DATAINloc,'RETURN_GAUSS_APPROX',0) ;



switch DATAINloc.GaussPointsToNodesMethod
    case 'GlobalSmoothinMassMatrix'
        for icomp = 1:DATAINloc.ncomponent
            VARNODES(icomp:ncomponent:end,:) = Nst\VARGAUSS(icomp:ncomponent:end,:)  ;
            if DATAINloc.RETURN_GAUSS_APPROX == 1
                VARGAUSS(icomp:ncomponent:end,:) = Nst*VARNODES(icomp:ncomponent:end,:)  ;
            end
        end
        
        
        
        
    otherwise
        error('Option  not implemented ')
end


