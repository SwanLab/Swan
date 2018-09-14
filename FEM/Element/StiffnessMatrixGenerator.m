classdef StiffnessMatrixGenerator < handle
    properties (Access = private)


        dimensionVariables
        
        
        dvolum
        
        ElementElastic
        Cmat

       
        StiffnesEntry
        StiffnesEntries
        IncrementalStifness
        Stifness
        
        changedElements
        ChangedEntriesIndices 
        StifnesEntriesChangedValues
        changedElementsRatio
        
        Bmat_all
        KeVector
        IncrementalStiffnessVector
        
        GlobalDofs
    end
    
    methods
        function obj= StiffnessMatrixGenerator(ElementElastic)
            obj.dimensionVariables = DimensionVariables();
            obj.storeDimensionalVariables(ElementElastic)
            
            obj.ElementElastic = ElementElastic;
            obj.initialize_dvolum(ElementElastic)

            connectivities = ElementElastic.geometry.interpolation.T;
            obj.GlobalDofs = GlobalDofs(obj.dimensionVariables,connectivities);
            obj.GlobalDofs.computeIJdofs()
            
            obj.initializeStifMat()
            obj.initializeStiffnesEntries()
            obj.storeBmat();
            obj.KeVector = StiffnessEntries(obj.dimensionVariables,obj.Bmat_all,obj.dvolum);
        end
        
        function generate(obj,Cmat)
             obj.obtainChangedElements(Cmat);
          
            if obj.changedElementsRatio ~= 0

            StiffnessIncrementEntries = obj.KeVector.compute(obj.Cmat,obj.changedElements);
            
            obj.computeIncrementalStiffnessVector(StiffnessIncrementEntries)
            obj.StiffnesEntries(obj.ChangedEntriesIndices) = StiffnessIncrementEntries;
            obj.assemble_matrix()
            
            obj.add_matrix()
            obj.symmetrizeStiffMat()
           end
        end
        
        function K = getStiffMatrix(obj)
            K = obj.Stifness;
        end
        
    end
    
    
    methods (Access = private)
    
        
        function computeIncrementalStiffnessVector(obj,StiffnesEntriesYY)
            obj.IncrementalStiffnessVector = StiffnesEntriesYY - obj.StiffnesEntries(obj.ChangedEntriesIndices);
        end
        
        function storeBmat(obj)
            obj.Bmat_all = zeros(obj.dimensionVariables.ngaus,obj.dimensionVariables.nstre,obj.dimensionVariables.nnode*obj.dimensionVariables.nunkn,obj.dimensionVariables.nelem);
            for igaus = 1:obj.dimensionVariables.ngaus
                obj.Bmat_all(igaus,:,:,:) = obj.computeBmatrix(igaus);
            end
        end
        
        
        function  obtainChangedElements(obj,Cmat)
            if ~isempty(obj.Cmat)
               
               dif = squeeze(Cmat(1,1,:) - obj.Cmat(1,1,:));
               for i= 1:size(Cmat,1)
                   for j = 1:size(Cmat,2)
                       dCmat = squeeze((Cmat(i,j,:) - obj.Cmat(i,j,:)));
                       dif = dif + dCmat.^2;
                   end
               end
               
               obj.changedElements = abs(dif)/norm(Cmat(:)) > 1e-15;
               obj.changedElementsRatio = sum(obj.changedElements)/size(obj.changedElements,1);
               obj.ChangedEntriesIndices = repmat(obj.changedElements,obj.dimensionVariables.ndofPerElement*obj.dimensionVariables.ndofPerElement,1);
              
            
            else
               obj.changedElements = true(obj.dimensionVariables.nelem,1);
               obj.ChangedEntriesIndices = repmat(obj.changedElements,obj.dimensionVariables.ndofPerElement*obj.dimensionVariables.ndofPerElement,1);
               obj.changedElementsRatio = 1;
               
            end
            obj.Cmat = Cmat;
               
        end
        
        function initializeStifMat(obj)
            obj.Stifness = sparse(obj.dimensionVariables.ndof,obj.dimensionVariables.ndof);
        end
        
        
        
        function storeDimensionalVariables(obj,ElementElastic)
            obj.dimensionVariables.nnode = ElementElastic.nnode;
            obj.dimensionVariables.nunkn = ElementElastic.dof.nunkn;
            obj.dimensionVariables.nstre = ElementElastic.nstre;
            obj.dimensionVariables.ndof  = ElementElastic.dof.ndof;
            obj.dimensionVariables.nelem = ElementElastic.nelem;
            obj.dimensionVariables.ndofPerElement = obj.dimensionVariables.nnode*obj.dimensionVariables.nunkn;
            
            ndofT = obj.dimensionVariables.ndofPerElement;
            nelem = obj.dimensionVariables.nelem;
            nStiffnesEntries = ndofT*ndofT*nelem;
            
            obj.dimensionVariables.nStiffnesEntries = nStiffnesEntries;

        end

        
        function initialize_dvolum(obj,ElementElastic)
            ElementElastic.quadrature.computeQuadrature('LINEAR');
            ElementElastic.interpolation_u.computeShapeDeriv(ElementElastic.quadrature.posgp)
            ElementElastic.geometry.computeGeometry(ElementElastic.quadrature,ElementElastic.interpolation_u);
            obj.dimensionVariables.ngaus = ElementElastic.quadrature.ngaus;
            obj.dvolum = ElementElastic.geometry.dvolu;
        end
        
        function Bmat = computeBmatrix(obj,igaus)
            Bmat = obj.ElementElastic.computeB(igaus);
        end        
 
        
        function initializeStiffnesEntries(obj)
            obj.StiffnesEntries = zeros(obj.dimensionVariables.nStiffnesEntries,1);
        end
        

        function assemble_matrix(obj)
            deltaK = obj.IncrementalStiffnessVector;
            IG = obj.GlobalDofs.Idofs(obj.ChangedEntriesIndices);
            JG = obj.GlobalDofs.Jdofs(obj.ChangedEntriesIndices);
            ndof = obj.dimensionVariables.ndof;
            obj.IncrementalStifness = sparse(IG,JG,deltaK,ndof,ndof);
        end
        
        function add_matrix(obj)
            deltaK = obj.IncrementalStifness ;
            obj.Stifness = obj.Stifness + deltaK ;
        end
       

        function symmetrizeStiffMat(obj)
            obj.Stifness = 1/2 * (obj.Stifness + obj.Stifness');
        end
        
    end
    

end

