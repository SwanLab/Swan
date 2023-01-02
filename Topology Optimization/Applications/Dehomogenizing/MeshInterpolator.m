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
            obj.createFunctions();
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

% 
            s.coord = [1 0; 0 1; 0 0; 1 1];
            s.connec = [3 1 2; 1 4 2];
            m = Mesh(s);

            obj.coarseMesh = m;
        end

        function mFine = remesh(obj,mC)
            %mC = obj.coarseMesh;
            mC.computeEdges();
            e = mC.edges;

            eMesh = mC.computeEdgeMesh();


            allCoords(:,1) = obj.computeFunctionInNodesAndEdges(mC,mC.coord(:,1));
            allCoords(:,2) = obj.computeFunctionInNodesAndEdges(mC,mC.coord(:,2));
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
        end
        
        function createFunctions(obj)
            f   = @(x) x(:,1).*x(:,2);
            mCoarse   = obj.coarseMesh;
            fCoarseC  = obj.createFunction(mCoarse,f(mCoarse.coord),'P1');
            mFine     = obj.remesh(mCoarse);        

            [fFineC] = obj.createContinousFunctions(fCoarseC,mCoarse,mFine);
            
          
            mCoarse   = obj.coarseMesh;
            mCoarseD = mCoarse.createDiscontinuousMesh();

            s.fValues = f(mCoarseD.computeBaricenter()');%reshape(chi,[],3); 
            s.connec  = mCoarse.connec; 
            s.type    = mCoarse.type;
            f = P0Function(s);

            s.mesh    = mCoarse;
            s.connec  = mCoarse.connec;
            p = Projector_toP1Discontinuous(s);
            fCoarseD = p.project(f);

            mFineD = obj.remesh(mCoarseD); 
            fFineD = obj.createDiscontinousFunctions(fCoarseD,mCoarseD,mFineD);
        end

        function fFine = createContinousFunctions(obj,fCoarse,m,mFine)
            fV    = fCoarse.fValues;
            fAll  = obj.computeFunctionInNodesAndEdges(m,fV);
            fFine = obj.createFunction(mFine,fAll,'P1');
        end


        function fFine = createDiscontinousFunctions(obj,fCoarse,m,mFine)
            f = squeeze(fCoarse.fValues);
            f = f(:);
            fAllD = obj.computeFunctionInNodesAndEdges(m,f);                   
            for i= 1:3
                fAll(1,i,:) = fAllD(mFine.connec(:,i));
            end
            fFine = obj.createFunction(mFine,fAll,'P1Disc');
        end        


        function f = computeFunctionInNodesAndEdges(obj,m,fNodes)
            fEdges   = obj.interpolateFunctionInEdges(m,fNodes);
            f = [fNodes;fEdges];
        end

        function fE = interpolateFunctionInEdges(obj,m,fNodes)
            me = m.computeEdgeMesh();  
            fe = obj.createFunction(me,fNodes,'P1'); 
            q = Quadrature.set(m.type);
            q.computeQuadrature('CONSTANT');
            xV = q.posgp;
            fE = squeeze(fe.evaluate(xV));  
        end        

        function f = createFunction(obj,m,fValues,fType)
            s.type    = m.type;
            s.connec  = m.connec;
            s.fValues = fValues;     
            s.functionType    = fType;
            f = FeFunction.create(s);
            f.plot(m)
        end
        
    end
    
end