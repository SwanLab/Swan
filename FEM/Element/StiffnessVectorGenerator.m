classdef StiffnessVectorGenerator < handle
   
    properties
       values 
    end
    
    
    properties (Access = private)
        Entries
        changedElem
        nChangedElements
        
        Kvalues
        nKvalues
        PartialKvalues
        
        dim
        
        TotalB
        Dvolum
        
        PartialB
        dvolum
        CtimesB
        
        Cmat
    end
    
    methods
        
        function obj = StiffnessVectorGenerator(dim,TotalB,TotalDvolum)
            obj.dim         = dim;
            obj.TotalB      = TotalB;
            obj.Dvolum      = TotalDvolum;
        end
        
        function compute(obj,DeltaC)
            obj.init(DeltaC)
            for igaus=1:obj.dim.ngaus
                obj.computeDvolum(igaus)
                obj.computePartialB(igaus)
                obj.computeCtimesBmatrix();
                obj.computeStiffnessEntries();
            end
            obj.values = obj.Kvalues;
        end
        
        function init(obj,DeltaC)
            obj.changedElem      = DeltaC.changedElements.values;
            obj.nChangedElements = DeltaC.changedElements.number;
            obj.Cmat             = DeltaC.Cmat(:,:,obj.changedElem);
            obj.Entries          = EntryIndex(obj.nChangedElements);
            
            obj.compute_nKvalues()
            obj.initializeKvalues();
        end
        
        function computeDvolum(obj,igaus)
            obj.dvolum =  obj.Dvolum(obj.changedElem,igaus);
        end
        
        function computePartialB(obj,igaus)
           obj.PartialB = squeeze(obj.TotalB(igaus,:,:,obj.changedElem));
        end
        
        function computeStiffnessEntries(obj)
            obj.Entries.init()
            for idof=1:obj.dim.ndofPerElement
                obj.computeDiagonalEntries(idof);
                for jdof=1:idof-1
                    obj.computeUpperAndLowerEntries(idof,jdof);
                end
            end
            obj.Kvalues;
        end
        
         function compute_nKvalues(obj)
            ndofs = obj.dim.ndofPerElement;
            obj.nKvalues = ndofs*ndofs*obj.nChangedElements;
        end
        
        
        function computeDiagonalEntries(obj,idof)
            obj.computeStiffEntries(idof,idof);
            obj.storeValueEntries();
        end
        
        function computeUpperAndLowerEntries(obj,idof,jdof)
            obj.computeStiffEntries(idof,jdof);
            obj.storeValueEntries();
            obj.storeValueEntries();
        end
        
        function storeValueEntries(obj)
            obj.Entries.compute();
            obj.addKvalues();
        end
        
        function computeCtimesBmatrix(obj)
            CB = zeros(obj.dim.nstre,obj.dim.ndofPerElement,obj.nChangedElements);
            for i=1:obj.dim.nstre
                PermutedCmat = permute(obj.Cmat(i,:,:),[2,1,3]);
                ReplicatedCmat = repmat(PermutedCmat,1,obj.dim.ndofPerElement,1);
                CB(i,:,:) = sum(ReplicatedCmat.* obj.PartialB,1);
            end
            obj.CtimesB = CB;
        end
        
        function initializeKvalues(obj)
            obj.Kvalues = zeros(obj.nKvalues,1);
        end
        
        function  computeStiffEntries(obj,idof,jdof)
            dV  = obj.dvolum;
            B   = obj.PartialB(:,idof,:);
            CB  = obj.CtimesB(:,jdof,:);
            BCB = sum( B.*CB ,1);
            BCB = squeeze(BCB);
            K   = BCB.*dV;
            obj.PartialKvalues = K;
        end
        
        function addKvalues(obj)
            In = obj.Entries.values;
            Kp = obj.PartialKvalues;
            Kt = obj.Kvalues;
            Kt(In,1) = Kt(In,1) + Kp;
            obj.Kvalues = Kt;
        end
        
    end
    
    
end

