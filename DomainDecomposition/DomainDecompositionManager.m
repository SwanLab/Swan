classdef DomainDecompositionManager < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        meshReference
        meshSubDomain
        nSubdomains
    end
    
    methods (Access = public)
        
        function obj = DomainDecompositionManager()
            obj.init();
            obj.createReferenceMesh()
            obj.createSubDomainMeshes();
            % obj.createSubDomainStiffnessMatrices();
            % obj.computeGlobalStiffnessMatrix();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.nSubdomains = [3 3];
        end
        
        function createReferenceMesh(obj)
            filename='lattice_ex1';
            a.fileName=filename;
            femD = FemDataContainer(a);   
            mS = femD.mesh;
            obj.meshReference = mS;
        end
        
        function L = computeReferenceMeshLength(obj) 
            coord = obj.meshReference.coord;
            Lx = max(coord(:,1));
            Ly = max(coord(:,2));            
            L = [Lx Ly];          
        end
        
        function createSubDomainMeshes(obj)
            nX = obj.nSubdomains(1);
            nY = obj.nSubdomains(2);   
            figure(2)
            for jDom = 1:nY
                for iDom = 1:nX
                    coordIJ = obj.computeSubdomainCoords(jDom,iDom);
                    mIJ     = obj.createSubdomainMesh(coordIJ);
                    mIJ.plot();
                    mD{iDom,jDom} = mIJ;
                    hold on
                end
            end
           obj.meshSubDomain = mD;
        end
        
        function m = createSubdomainMesh(obj,coord)
            connec0  = obj.meshReference.connec;            
            s.coord  = coord;
            s.connec = connec0;
            m = Mesh(s);
        end
        
        function coord = computeSubdomainCoords(obj,jDom,iDom)
            coord0  = obj.meshReference.coord;  
            L = obj.computeReferenceMeshLength();
            Lx = L(1);
            Ly = L(2);
            coord(:,1) = coord0(:,1)+Lx*(iDom-1);
            coord(:,2) = coord0(:,2)+Ly*(jDom-1);
        end
        
        function gConnec = computeGlobalConnec(obj,connec)
            nnodes    = obj.meshReference.nnodes;
            nX = obj.nSubdomains(1);            
            gConnec = nnodes*(nX*(jDom-1)+iDom-1);            
        end
        
    end
end
        