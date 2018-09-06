classdef StiffnessMatrixGenerator < handle 
    properties (Access = private)
        nnode
        nunkn
        nstre
        ndof
        nelem
        ndofPerElement
        inodes
        icomps
        ngaus
        
        connectivities
        dvolum

        ElementElastic
        Cmat
        Bmat
        CtimesB
        
        Global_Idof
        Global_Jdof
        Global_Idofs
        Global_Jdofs
        NumberOfGlobalDofs
        InitialEntryIndex
        FinalEntryIndex
        
        StiffnesEntry
        StiffnesEntries        
        Subintegrated_StifMat
        StifMat
    end
    
    methods
        function obj= StiffnessMatrixGenerator(ElementElastic)
            obj.storeDimensionalVariables(ElementElastic)
            obj.connectivities = ElementElastic.geometry.interpolation.T;
            obj.Cmat = ElementElastic.material.C;
            obj.ElementElastic = ElementElastic;
            obj.initialize_dvolum(ElementElastic)
        end
        
        function generate(obj)
            obj.initializeStifMat()
            
            for igaus=1:obj.ngaus
                obj.initializeGlobal_IJdofs();
                obj.computeBmatrix(igaus);
                obj.compute_CtimesBmatrix();
                obj.computeMatrixEntries(igaus);
                obj.assemble_matrix()
                obj.add_matrix()
            end
            obj.symmetrizeStiffMat()
        end
        
        function K = getStiffMatrix(obj)
            K = obj.StifMat;
        end
        
    end
    
    methods (Access = private)
        function initializeStifMat(obj)
            obj.StifMat = sparse(obj.ndof,obj.ndof);
        end
        
        function storeDimensionalVariables(obj,ElementElastic)
            obj.nnode = ElementElastic.nnode;
            obj.nunkn = ElementElastic.dof.nunkn;
            obj.nstre = ElementElastic.nstre;
            obj.ndof  = ElementElastic.dof.ndof;
            obj.nelem = ElementElastic.nelem;
            obj.ndofPerElement = obj.nnode*obj.nunkn;
            obj.inodes=reshape(repmat(1:obj.nnode,obj.nunkn,1),1,[]);
            obj.icomps=repmat(1:obj.nunkn,1,obj.nnode);
        end
        
        
        function initialize_dvolum(obj,ElementElastic)
            ElementElastic.quadrature.computeQuadrature('LINEAR');
            ElementElastic.interpolation_u.computeShapeDeriv(ElementElastic.quadrature.posgp)
            ElementElastic.geometry.computeGeometry(ElementElastic.quadrature,ElementElastic.interpolation_u);
            obj.ngaus = ElementElastic.quadrature.ngaus;
            obj.dvolum = ElementElastic.geometry.dvolu;
        end
        
        function computeBmatrix(obj,igaus)
            obj.Bmat = obj.ElementElastic.computeB(igaus);
        end
        

        
        function initializeGlobal_IJdofs(obj)
            obj.Global_Idofs = obj.initializeGlobalDofs();
            obj.Global_Jdofs = obj.initializeGlobalDofs();
        end
        
        function GlobalDofs = initializeGlobalDofs(obj)
            GlobalDofs = zeros(obj.ndofPerElement*obj.ndofPerElement*obj.nelem,1);
        end
        
        function computeMatrixEntries(obj,igaus)

            obj.initializeGaussLoopVariables()
            for idof=1:obj.ndofPerElement
                obj.obtainGlobal_Idof(idof);
                obj.obtainNumberOfGlobalDofs();
                obj.computeDiagonalEntries(idof,igaus)
                for jdof=1:idof-1
                    obj.obtainGlobal_Jdof(jdof);
                    obj.computeUpperAndLowerEntries(idof,jdof,igaus);
                end
            end
        end
        
        function initializeGaussLoopVariables(obj)
            obj.initializeInitialEntryIndex()
            obj.initializeStiffnesEntries()
        end
        
        function computeDiagonalEntries(obj,idof,igaus)
            obj.computeStiffEntries(idof,idof,igaus);
            obj.storeDiagonalValues()
        end
        
        function computeUpperAndLowerEntries(obj,idof,jdof,igaus)
            obj.computeStiffEntries(idof,jdof,igaus);
            obj.storeUpperDiagonalValues()
            obj.storeLowerDiagonalValues()
        end
        
        
        function storeDiagonalValues(obj)
            obj.storePositionAndValueEntries(obj.Global_Idof,obj.Global_Idof)
        end
        
        function storeUpperDiagonalValues(obj)
            obj.storePositionAndValueEntries(obj.Global_Idof,obj.Global_Jdof)
        end
        
        function storeLowerDiagonalValues(obj)
            obj.storePositionAndValueEntries(obj.Global_Jdof,obj.Global_Idof)
        end
        
        function storePositionAndValueEntries(obj,Global_Idof,Global_Jdof)
            obj.updateFinalEntryIndex()
            EntriesIndex = obj.obtainEntriesIndex();
            obj.Global_Idofs(EntriesIndex,1)    =  Global_Idof;
            obj.Global_Jdofs(EntriesIndex,1)    =  Global_Jdof;
            obj.StiffnesEntries(EntriesIndex,1) =  obj.StiffnesEntry;
            obj.updateInitialEntryIndex()
        end
        
        function initializeInitialEntryIndex(obj)
            obj.InitialEntryIndex=1;
        end
        
        function updateFinalEntryIndex(obj)
            obj.FinalEntryIndex = obj.InitialEntryIndex + obj.NumberOfGlobalDofs -1;
        end
        
        function updateInitialEntryIndex(obj)
            obj.InitialEntryIndex = obj.FinalEntryIndex + 1;
        end
        
        function EntriesIndex = obtainEntriesIndex(obj)
            EntriesIndex = obj.InitialEntryIndex:obj.FinalEntryIndex;
        end
                
        function obtainGlobal_Idof(obj,idof)
            obj.Global_Idof = obj.transformLocal2Global(idof);
        end
        
        function obtainGlobal_Jdof(obj,jdof)
            obj.Global_Jdof = obj.transformLocal2Global(jdof);
        end
                
        function GlobalDofs = transformLocal2Global(obj,LocalDof)
            GlobalDofs = obj.nunkn*(obj.connectivities(:,obj.inodes(LocalDof))-1)+obj.icomps(LocalDof);
        end
      
        function obtainNumberOfGlobalDofs(obj)
            obj.NumberOfGlobalDofs = length(obj.Global_Idof);
        end
        
        function compute_CtimesBmatrix(obj)
            CB = zeros(obj.nstre,obj.ndofPerElement,obj.nelem);
            for i=1:obj.nstre
                PermutedCmat = permute(obj.Cmat(i,:,:),[2,1,3]);
                ReplicatedCmat = repmat(PermutedCmat,1,obj.ndofPerElement,1);
                CB(i,:,:) = sum(ReplicatedCmat.* obj.Bmat,1);
            end
            obj.CtimesB = CB;
        end

        function initializeStiffnesEntries(obj)
            obj.StiffnesEntries = zeros(obj.ndofPerElement*obj.ndofPerElement*obj.nelem,1);
        end       
        
        function computeStiffEntries(obj,idof,jdof,igaus)
            BCB = squeeze(sum(obj.Bmat(:,idof,:) .* obj.CtimesB(:,jdof,:),1));
            obj.StiffnesEntry = obj.dvolum(:,igaus).*BCB;
        end
        
        function assemble_matrix(obj)
            obj.Subintegrated_StifMat = sparse(obj.Global_Idofs,obj.Global_Jdofs,obj.StiffnesEntries,obj.ndof,obj.ndof);
        end
        
        function add_matrix(obj)
            obj.StifMat = obj.StifMat + obj.Subintegrated_StifMat;
        end
        
        function symmetrizeStiffMat(obj)
            obj.StifMat = 1/2 * (obj.StifMat + obj.StifMat');
        end  
    end
end

