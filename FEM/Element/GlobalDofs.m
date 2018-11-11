classdef GlobalDofs < handle
    
    properties
        Idofs
        Jdofs
    end
    
    properties (Access = private)
        Idof
        Jdof
        
        nodesInElement
        VectorDimensions
        nunkn
        dofsPerElement
        connectivities
        
        EntryIndex
    end
    
    methods (Access = public)
        
        function obj = GlobalDofs(connectivities,dim)
            obj.connectivities   = connectivities;
            obj.EntryIndex       = EntryIndex(dim.nelem);
            obj.nunkn            = dim.nunkn;
            obj.dofsPerElement   = dim.ndofPerElement;
            obj.nodesInElement   = reshape(repmat(1:dim.nnode,dim.nunkn,1),1,[]);
            obj.VectorDimensions = repmat(1:dim.nunkn,1,dim.nnode);
        end
        
       
        function computeIJdofs(obj)
            obj.EntryIndex.init()
            for idof=1:obj.dofsPerElement
                obj.obtainIdof(idof);
                obj.storeDiagonalEntries()
                for jdof=1:idof-1
                    obj.obtainJdof(jdof);
                    obj.storeUpperDiagonalEntries()
                    obj.storeLowerDiagonalEntries()
                end
            end            
        end
        
    end
    
    methods (Access = private)
        
        function storeDiagonalEntries(obj)
            obj.storePositionEntries(obj.Idof,obj.Idof)
        end
        
        function storeUpperDiagonalEntries(obj)
            obj.storePositionEntries(obj.Idof,obj.Jdof)
        end
        
        function storeLowerDiagonalEntries(obj)
            obj.storePositionEntries(obj.Jdof,obj.Idof)
        end
        
        function storePositionEntries(obj,Idof,Jdof)
            obj.EntryIndex.compute();
            obj.Idofs(obj.EntryIndex.values,1) =  Idof;
            obj.Jdofs(obj.EntryIndex.values,1) =  Jdof;
        end
        
        function obtainIdof(obj,idof)
            obj.Idof = obj.transformLocal2Global(idof);
        end
        
        function obtainJdof(obj,jdof)
            obj.Jdof = obj.transformLocal2Global(jdof);
        end
        
        function GlobalDofs = transformLocal2Global(obj,LocalDof)
            LocalNode       = obj.nodesInElement(LocalDof);
            VectorDimension = obj.VectorDimensions(LocalDof);
            GlobalDofs = obj.nunkn*(obj.connectivities(:,LocalNode)-1) + VectorDimension;
        end 
        
        
    end
    
    
end

