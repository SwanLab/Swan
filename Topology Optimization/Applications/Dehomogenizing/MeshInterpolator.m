classdef MeshInterpolator < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
       fineMesh 
       coarseMesh
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = MeshInterpolator()
            close all
            obj.createCoarseMesh()
            obj.createFineMesh();
            obj.createCoarseFeFunction();
            obj.interpolateFeFunction();
        end
        
    end
    
    methods (Access = private)
        
        
        function createCoarseMesh(obj)
            s.xmin = 0;
            s.xmax = 1;
            s.ymin = 0;
            s.ymax = 1;
            s.nInteriorPoints = 500;
            s.nBoundaryPoints = 15;
            uM = UnstructuredMeshCreator(s);
            m = uM.create();
            obj.coarseMesh = m;
        end
        
        function createFineMesh(obj)
            mC = obj.coarseMesh;
            mC.computeEdges();
            e = mC.edges;

            eMesh = mC.computeEdgeMesh();

            newCoord = eMesh.computeBaricenter()';
            allCoord = [mC.coord;newCoord];

            
            nnodes            = eMesh.nnodes;
            newNodes(:,1)     = (nnodes+1):(nnodes+eMesh.nelem);
            edgesWithNewNodes = 1:eMesh.nelem;
            edgesAsBefore     = setdiff(1:eMesh.nelem,edgesWithNewNodes);
            firstElems    = [eMesh.connec(edgesWithNewNodes,1),newNodes];
            secondElems   = [newNodes,eMesh.connec(edgesWithNewNodes,2)];
            newAllElems   = [firstElems,secondElems];
            newConecs     = reshape(newAllElems(:)',[],2);

            allConnecEdges = [eMesh.connec(edgesAsBefore,:),newConecs];                        

            sE.coord  = allCoord;
            sE.connec = allConnecEdges;
            sE.kFace  = eMesh.kFace;
            eMeshFine = Mesh(sE);
            figure()
            eMeshFine.plot()



            rCells = 1:mC.nelem;
            fCells = setdiff(1:mC.nelem,rCells);

            nV1 = mC.connec(rCells,1);
            nV2 = mC.connec(rCells,2);
            nV3 = mC.connec(rCells,3);
            
            nE1 = newNodes(e.edgesInElem(rCells,1),1);
            nE2 = newNodes(e.edgesInElem(rCells,2),1);
            nE3 = newNodes(e.edgesInElem(rCells,3),1);

            nConnec(:,1,:) = [nV1 nE1 nE3];
            nConnec(:,2,:) = [nE1 nV2 nE2];
            nConnec(:,3,:) = [nE1 nE2 nE3];
            nConnec(:,4,:) = [nE3 nE2 nV3];

            newConnec = reshape(nConnec,[],3);
            
            allConnec = [mC.connec(fCells,:),newConnec]; 

            figure()
            sF.connec = allConnec;
            sF.coord  = [mC.coord;newCoord];
            mFine = Mesh(sF);
            mFine.plot()
            obj.fineMesh = mFine;
        end

        function createCoarseFeFunction(obj)
            f = @(x) x(:,1).*x(:,2);
            m = obj.coarseMesh;
            s.connec  = m.connec;
            s.type    = m.type;
            s.fValues = f(m.coord);
            cF = P1Function(s); 
            cF.plot(m)
        end

        function interpolateFeFunction(obj)
            f = @(x) x(:,1).*x(:,2);            
            m = obj.coarseMesh.computeEdgeMesh();  
            fValues = f(m.coord);
            s.connec  = m.connec;
            s.type    = m.type;
            s.fValues = fValues;
            cF = P1Function(s);             
            q = Quadrature.set(m.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            fNewValues = squeeze(cF.evaluate(xV));                


            allF = [fValues;fNewValues];


            m = obj.fineMesh;
            s.connec  = m.connec;
            s.type    = m.type;
            s.fValues = allF;
            cF = P1Function(s); 
            cF.plot(m)            

        end
        
    end
    
end