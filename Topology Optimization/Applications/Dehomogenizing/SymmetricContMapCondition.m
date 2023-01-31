classdef SymmetricContMapCondition < handle
    
    properties (Access = private)
        isCoherent
        orientationDisc        
    end
    
    properties (Access = private)
       meshDisc
       meshCont
       orientation
    end
    
    methods (Access = public)
        
        function obj = SymmetricContMapCondition(cParams)
            obj.init(cParams);            
        end
        
        function c = computeCondition(obj)
            obj.createOrientationDiscontinous();
            obj.isOrientationCoherent();
            c = obj.computeSymmetricCondition();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)         
            obj.meshCont    = cParams.meshCont;
            obj.meshDisc    = cParams.meshDisc;
            obj.orientation = cParams.orientation;
        end
        
        function createOrientationDiscontinous(obj)
            s.connec = obj.meshCont.connec;
            s.type   = obj.meshCont.type;
            s.fNodes = obj.orientation;
            
            s.fValues = obj.orientation;
            s.connec = obj.meshCont.connec;
            s.type   = obj.meshCont.type;
            f = P1Function(s);
            s.mesh   = obj.meshCont;
            s.connec = obj.meshCont.connec;
            p = Projector_toP1Discontinuous(s);
            fD = p.project(f);
            
            obj.orientationDisc = fD.fValues;
        end         

        function isOrientationCoherent(obj)
            s.mesh        = obj.meshDisc;
            s.orientation = obj.orientationDisc;
            c = CoherentOrientationSelector(s);
            isC = c.isOrientationCoherent();
            obj.isCoherent = isC;
        end

        function sC = computeSymmetricCondition(obj)
            nnodeD    = obj.meshDisc.nnodeElem;
            nElemD    = obj.meshDisc.nelem;                                    
            nnodesC   = obj.meshCont.nnodes;    
            connecC   = obj.meshCont.connec;
            connecD   = obj.meshDisc.connec;
            sC = sparse(nnodeD*nElemD,nnodesC);
            for iNode = 1:nnodeD
                isC  = obj.isCoherent(:,iNode);
                cond = obj.computeConformalMapCondition(isC);
                nodesC = connecC(:,iNode);
                nodesD = connecD(:,iNode);
                sC = sC + sparse(nodesD,nodesC,cond,nnodeD*nElemD,nnodesC);
            end
        end
       
    end

    methods (Access = private, Static)

        function m = computeConformalMapCondition(isCoherent)
            isMapSymmetric = isCoherent;
            isMapAntiSymm  = ~isCoherent;
            m = zeros(size(isMapSymmetric));
            m(isMapSymmetric) = 1;
            m(isMapAntiSymm)  = -1;
        end        

    end
    
end