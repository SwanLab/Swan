classdef StiffnessMatrixGenerator < handle
    properties (Access = private)
        
        ndof
        nentries
        
        Cmat
        DeltaC
        
        K
        DeltaK
        
        Kvector
        DeltaKvector

        
        KvectorGenerator
        
        DofsGlobal
    end
    
    methods
        function obj= StiffnessMatrixGenerator(connec,Bmat,dvolum,dim)
            obj.ndof     = dim.ndof;
            obj.nentries = dim.nentries;
            
            obj.DeltaC = ConstitutiveTensorIncrement(dim);
            obj.computeGlobalDofs(connec,dim)

            obj.initK()
            obj.initKvector()
            
            obj.KvectorGenerator = StiffnessVectorGenerator(dim,Bmat,dvolum);
        end
        
        function generate(obj,Cmat)
            
            obj.DeltaC.obtainChangedElements(Cmat)
            
            if obj.DeltaC.changedElements.Ratio ~= 0
                
                obj.KvectorGenerator.compute(obj.DeltaC);
               
                obj.computeDeltaKvector()
                obj.updateKvector()
                obj.assemble_matrix()
                
                obj.add_matrix()
                obj.symmetrizeStiffMat()
            end
        end
        
        function K = getStiffMatrix(obj)
            K = obj.K;
        end
        
    end
    
    
    methods (Access = private)
        
        function  updateKvector(obj)
            Index = obj.DeltaC.changedElements.Entries;
            obj.Kvector(Index) = obj.KvectorGenerator.values;
        end
        
        
        function computeGlobalDofs(obj,connec,dim)
            obj.DofsGlobal = GlobalDofs(connec,dim);
            obj.DofsGlobal.computeIJdofs()
        end
        
        function computeDeltaKvector(obj)
            obj.DeltaKvector = obj.KvectorGenerator.values - obj.Kvector(obj.DeltaC.changedElements.Entries);
        end
        
        function initK(obj)
            obj.K = sparse(obj.ndof,obj.ndof);
        end
               
        function initKvector(obj)
            obj.Kvector = zeros(obj.nentries,1);
        end
        
        function assemble_matrix(obj)
            IG = obj.DofsGlobal.Idofs(obj.DeltaC.changedElements.Entries);
            JG = obj.DofsGlobal.Jdofs(obj.DeltaC.changedElements.Entries);
            obj.DeltaK = sparse(IG,JG,obj.DeltaKvector,obj.ndof,obj.ndof);
        end
        
        function add_matrix(obj)
            obj.K = obj.K + obj.DeltaK ;
        end
        
        function symmetrizeStiffMat(obj)
            obj.K = 1/2 * (obj.K + obj.K');
        end
        
    end
    
    
end

