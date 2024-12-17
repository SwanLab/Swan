classdef OrientedMappingComputer < handle

    properties (Access = private)  
        orientation        
        isCoherent        
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
        
        function obj = OrientedMappingComputer(cParams)
            obj.init(cParams)
        end

        function dCoord = computeDeformedCoordinates(obj)
            obj.dilation              = obj.computeDilation();
            obj.dilatedOrientation    = obj.computeDilatedOrientationVector();            
            b1 = obj.orientation{1};
            obj.isCoherent            = obj.computeIsOrientationCoherent(b1);
            obj.interpolator          = obj.computeInterpolator();
            obj.phiMapping            = obj.computeMappings();
            obj.totalCorrector        = obj.computeTotalCorrector(b1);
            dCoord                    = obj.phiMapping + obj.totalCorrector;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.orientation = cParams.orientation; %%importnat since regul need P1
        end

        function isCF = computeIsOrientationCoherent(obj,b1)

            % s.filterType = 'LUMP';
            % s.mesh  = obj.mesh;
            % s.trial = LagrangianFunction.create(obj.mesh,a1.ndimf,'P1');
            % f = Filter.create(s);
            % a1 = f.compute(a1,2);


      %      a1   = project(a1,'P1'); 
       %     u = UnitBallProjector([]);
            
      %      a1.setFValues(u.project(a1.fValues))

            b1   = b1.getFvaluesByElem();
            bN1   = squeeze(b1(:,1,:));            
            isCoh = false(obj.mesh.nnodeElem,obj.mesh.nelem);
            for iNode = 1:obj.mesh.nnodeElem
                bNi    = squeeze(b1(:,iNode,:));
                bN1bNI = dot(bN1,bNi);
                isCoh(iNode,:) = (bN1bNI)>0;
            end
            s.fValues = isCoh(:);
            s.mesh    = obj.mesh;
            s.order   = 'P1D';
            isCF = LagrangianFunction(s); 
        end

        function sC = computeInterpolator(obj)
            nnodeD    = obj.mesh.nnodeElem;
            nElemD    = obj.mesh.nelem;
            nnodesC   = obj.mesh.nnodes;
            connecC   = obj.mesh.connec;
            connecD   = obj.isCoherent.getDofConnec();
            sC = sparse(nnodeD*nElemD,nnodesC);
            isC = obj.isCoherent.getFvaluesByElem();
            for iNode = 1:nnodeD
                isCi   = squeeze(isC(1,iNode,:));
                cond   = obj.computeConformalMapCondition(isCi);
                nodesC = connecC(:,iNode);
                nodesD = connecD(:,iNode);
                sC = sC + sparse(nodesD,nodesC,cond,nnodeD*nElemD,nnodesC);
            end
        end

        function isSf = computeSingularities(obj,a1)
            a1 = project(a1,'P1D');
            aD = permute(a1.getFvaluesByElem(), [1 3 2]);
            a1 = zeros(3,obj.mesh.nelem);
            a2 = zeros(3,obj.mesh.nelem);
            a3 = zeros(3,obj.mesh.nelem);
            a1(1:2,:) = aD(:,:,1);
            a2(1:2,:) = aD(:,:,2);
            a3(1:2,:) = aD(:,:,3);
            a1a2 = dot(a1,a2);
            a1a3 = dot(a1,a3);
            a2a3 = dot(a2,a3);
            isS = sign(a1a2.*a1a3.*a2a3)';
            s.fValues = isS<0;
            s.order   = 'P0';
            s.mesh    = obj.mesh;
            isSf = LagrangianFunction(s);
        end          

        function d = computeDilation(obj)
            s.orientationVector = obj.orientation;
            s.mesh              = obj.mesh;
            dC = DilationComputer(s);
            d  = dC.compute();
        end    

        function d = computeDilatedOrientationVector(obj)
            er = exp(obj.dilation);
            d = cell(obj.mesh.ndim,1);
            for iDim = 1:obj.mesh.ndim
                b  = obj.orientation{iDim};
                d{iDim} = er.*b;
            end
        end            

        function phiM = computeMappings(obj)
            s.mesh               = obj.mesh;
            s.interpolator       = obj.interpolator;
            s.dilatedOrientation = obj.dilatedOrientation;
            mC   = MappingComputer(s);
            phiM = mC.compute();
        end

        function tC = computeTotalCorrector(obj,a1)
           s.mesh = obj.mesh;
           s.isCoherent            = obj.isCoherent;
           s.isOrientationSingular = obj.computeSingularities(a1);
           s.interpolator          = obj.interpolator;
           s.dilatedOrientation    = obj.dilatedOrientation;
           s.phiMapping            = obj.phiMapping;
           tCC = TotalCorrectorComputer(s);           
           tC = tCC.compute();
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