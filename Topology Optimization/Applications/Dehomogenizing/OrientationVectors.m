classdef OrientationVectors < handle

    properties (Access = private)  
        value        
        isCoherent        
        singularities
        interpolator    
        dilation        
        dilatedOrientation               
        phiMapping
        totalCorrector
        phi        
    end
    
    properties (Access = private)
      mesh
      theta
    end
    
    methods (Access = public)
        
        function obj = OrientationVectors(cParams)
            obj.init(cParams)
        end

        function dCoord = computeDeformedCoordinates(obj)
            obj.createOrientationVector();
            obj.computeIsOrientationCoherent();
            obj.computeInterpolator();
            obj.computeSingularities();
            obj.computeDilation();
            obj.computeDilatedOrientationVector();
            obj.computeMappings();
            obj.computeTotalCorrector();
            obj.computeMappingWithSingularities();
            dCoord = obj.phi;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh  = cParams.mesh;
            obj.theta = cParams.theta;
        end
        
        function createOrientationVector(obj)
            a1(:,1) = cos(obj.theta);
            a1(:,2) = sin(obj.theta);
            a2(:,1) = -sin(obj.theta);
            a2(:,2) = cos(obj.theta);
            a(:,:,1) = a1;
            a(:,:,2) = a2;
            nDim = obj.mesh.ndim;
            orientation = cell(nDim,1);
            for iDim = 1:nDim
                s.fValues = a(:,:,iDim);
                s.mesh   = obj.mesh;
                bf = P1Function(s);
                orientation{iDim} = bf;
            end 
            obj.value = orientation;
        end
        
        function computeIsOrientationCoherent(obj)
            nnode = obj.mesh.nnodeElem;
            nElem = obj.mesh.nelem;
            isCoh = false(1,nnode,nElem);
            a1D   = obj.value{1}.project('P1D'); 
            a1    = a1D.fValues;
            aN1   = squeeze(a1(:,1,:));
            for iNode = 1:nnode
                aNi    = squeeze(a1(:,iNode,:));
                aN1aNI = dot(aN1,aNi);
                isCoh(1,iNode,:) = (aN1aNI)>0;
            end
            s.fValues = isCoh;
            s.mesh    = obj.mesh;
            isCF       = P1DiscontinuousFunction(s);            
            obj.isCoherent = isCF;
        end

        function sC = computeInterpolator(obj)
            nnodeD    = obj.mesh.nnodeElem;
            nElemD    = obj.mesh.nelem;
            nnodesC   = obj.mesh.nnodes;
            connecC   = obj.mesh.connec;
            connecD   = obj.isCoherent.computeDiscontinuousConnectivities();
            isCo = obj.isCoherent;
            sC = sparse(nnodeD*nElemD,nnodesC);
            for iNode = 1:nnodeD
                isC  = squeeze(isCo.fValues(1,iNode,:));
                cond = obj.computeConformalMapCondition(isC);
                nodesC = connecC(:,iNode);
                nodesD = connecD(:,iNode);
                sC = sC + sparse(nodesD,nodesC,cond,nnodeD*nElemD,nnodesC);
            end
            obj.interpolator = sC;
        end

        function computeSingularities(obj)
            s.mesh        = obj.mesh;
            s.orientation = obj.value{1};
            sC = SingularitiesComputer(s);
            sC.compute();
            obj.singularities = sC;%sCoord;
        end          

        function computeDilation(obj)
            s.orientationVector = obj.value;
            s.mesh              = obj.mesh;
            dC = DilationComputer(s);
            d  = dC.compute();
            obj.dilation = d;
        end    

        function computeDilatedOrientationVector(obj)
            s.fValues = exp(obj.dilation.fValues);
            s.mesh    = obj.mesh;
            er = P1Function(s);
            for iDim = 1:obj.mesh.ndim
                b  = obj.value{iDim};
                dO = P1Function.times(er,b);
                obj.dilatedOrientation{iDim} = dO;
            end
        end            

        function computeMappings(obj)
            s.mesh               = obj.mesh;
            s.interpolator       = obj.interpolator;
            s.dilatedOrientation = obj.dilatedOrientation;
            mC   = MappingComputer(s);
            phiM = mC.compute();
            obj.phiMapping = phiM;
        end

        function computeTotalCorrector(obj)
           s.mesh = obj.mesh;
           s.isCoherent         = obj.isCoherent;
           s.singularities      = obj.singularities;
           s.interpolator       = obj.interpolator;
           s.dilatedOrientation = obj.dilatedOrientation;
           s.phiMapping   = obj.phiMapping;
           tC = TotalCorrectorComputer(s);           
           obj.totalCorrector = tC.compute();
        end

        function computeMappingWithSingularities(obj)
            psiTs = P1DiscontinuousFunction.sum(obj.phiMapping,obj.totalCorrector);
            obj.phi = psiTs;
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