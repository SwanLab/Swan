classdef GradientVariationWithBoundaryComputer < handle
    
    properties (Access = public)
       error 
       mean
       desv
    end
    
    properties (Access = private)
        levelSet
        curvature
        rPerimeter
        backgroundMesh
        boundaryMesh
        domainLength
        circleCase
        filePlotName        
        circleMesh
        unfittedMesh        
        integrator
        
        gradientCircunf
        nEpsilon
    end
    
    methods (Access = public)
        
        function obj = GradientVariationWithBoundaryComputer(cParams)
            obj.init(cParams)
        end
            
        function compute(obj)      
            obj.createIntegrator();
            for iEpsilon = 1:obj.nEpsilon
                obj.computeGradientInCircunference(iEpsilon);
                obj.computeVariables(iEpsilon);                
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.filePlotName   = cParams.filePlotName;
            obj.backgroundMesh = cParams.backgroundMesh;
            obj.circleMesh     = cParams.circleMesh;
            obj.levelSet       = cParams.levelSet;
            obj.curvature      = cParams.curvature;
            obj.rPerimeter     = cParams.regularizedPerimeter;
            obj.domainLength   = cParams.domainLength;
            obj.circleCase     = cParams.circleCase;
            obj.nEpsilon       = size(obj.rPerimeter.epsilons,2);
        end
                  
        function computeVariables(obj,iepsilon)
           int = obj.integrator;
           M = int.computeLHS();           
           k = obj.curvature;
           g = obj.gradientCircunf;
           obj.error(iepsilon) = (g-k)'*M*(g-k)/(k*sum(M(:))*k);
           obj.mean(iepsilon) = sum(int.integrate(g))/sum(int.integrate(ones(size(g))));
           obj.desv(iepsilon) = (g-obj.mean(iepsilon))'*M*(g-obj.mean(iepsilon))/(sum(M(:)));
        end
        
        function createIntegrator(obj)
           m = obj.circleMesh.mesh;
           m = m.computeCanonicalMesh(); 
           s.type = 'SIMPLE';
           s.mesh = m;           
           s.npnod = m.nnodes;
           s.globalConnec = m.connec;
           obj.integrator = Integrator.create(s);            
        end        
        
        function computeGradientInCircunference(obj,iepsilon)
            nCell = obj.circleMesh.mesh.nelem;
            gB = zeros(nCell,1);                                    
            gD = obj.rPerimeter.perimetersGradient(:,iepsilon);
            switch obj.circleCase
                case 'interior'
                    for icell = 1:nCell
                        gNodalInCell  = obj.computeFnodalInCell(icell,gD);
                        xPos          = obj.computeXpos(icell);
                        gInXpos       = obj.interpolateValue(xPos,gNodalInCell);
                        gB(icell) = gInXpos;
                    end
                case 'exterior'
                    m = obj.circleMesh;
                    gB = gD(m.nodesInBoxFaces,1);
            end
            obj.gradientCircunf = gB;
        end        
        
        function fNodalInCell = computeFnodalInCell(obj,icell,fNodal)
            m = obj.circleMesh;            
            cellGlobal   = m.cellContainingSubcell(icell);
            nodeTriangle = obj.backgroundMesh.connec(cellGlobal,:);
            fNodalInCell = fNodal(nodeTriangle(:),:);
        end
        
        function x1pos = computeXpos(obj,icell)
            m = obj.circleMesh;
            x1pos = m.xCoordsIso(:,1,icell);
            x1pos = squeeze(x1pos);
        end
        
        function fInterp = interpolateValue(obj,xpos,fNodal)
            interp = Interpolation.create(obj.backgroundMesh.type,'LINEAR');
            shape = interp.computeShapeFunctions(xpos);
            fInterp = 0;
            for inode = 1:length(shape)
                fInterp = fInterp + fNodal(inode)*shape(inode);
            end
        end
        
    end    
  
end