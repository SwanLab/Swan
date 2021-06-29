classdef KGeneratorWithfullStoredB < handle
    
    properties
        nodesInElement
        VectorDimensions
        dofsPerElement
        connectivities
        Btot
        Bfull
        nunkn
        dim
        ndofGlobal
        CmatTot
        K
        dvolum
        nt
    end
    
    methods
        
        function obj = KGeneratorWithfullStoredB(dim,connectivities,Bfull,dvolum)
            obj.dofsPerElement   = dim.ndofPerElement;
            obj.nodesInElement   = reshape(repmat(1:dim.nnode,dim.nunkn,1),1,[]);
            obj.VectorDimensions = repmat(1:dim.nunkn,1,dim.nnode);
            obj.connectivities   = connectivities;
            obj.Bfull = Bfull;
            obj.dim = dim; 
            obj.nunkn = dim.nunkn;
            obj.dvolum = dvolum;
            
            obj.ndofGlobal = max(max(connectivities))*dim.nunkn;
                
            obj.nt = obj.dim.ngaus*obj.dim.nelem*obj.dim.nstre;    
            obj.computeBtot();  
        end
        
        function compute(obj,Cmat)

                      
            obj.computeCmatBlockDiagonal(Cmat);
            obj.computeStiffnes()            
        end
        
        function  computeStiffnes(obj)

            B = obj.Btot;
            CB = obj.CmatTot*B;
            obj.K = B'*CB;            
        end
        
        
        function computeCmatBlockDiagonal(obj,Cmat)
            obj.CmatTot = sparse(obj.nt,obj.nt);
            for istre = 1:obj.dim.nstre
                for jstre = 1:obj.dim.nstre
                    for igaus = 1:obj.dim.ngaus
                        posI = (istre)+(obj.dim.nstre)*(igaus-1) : obj.dim.ngaus*obj.dim.nstre : obj.nt ;
                        posJ = (jstre)+(obj.dim.nstre)*(igaus-1) : obj.dim.ngaus*obj.dim.nstre : obj.nt ;
                        
                        Ct = squeeze(Cmat(istre,jstre,:,igaus)).*obj.dvolum(:,igaus);                        
                        obj.CmatTot = obj.CmatTot + sparse(posI,posJ,Ct,obj.nt,obj.nt);
                    end
                        
                end
            end
        end
        
        function GlobalDofs = transformLocal2Global(obj,LocalDof)
            LocalNode        = obj.nodesInElement(LocalDof);
            VectorDimension  = obj.VectorDimensions(LocalDof);
            GlobalDofs       = obj.nunkn*(obj.connectivities(:,LocalNode)-1) + VectorDimension;
        end
        
        
        
        function GlobalDofs = computeBtot(obj)
            obj.Btot = sparse(obj.nt,obj.ndofGlobal);
            for idof=1:obj.dofsPerElement
                GlobalDofs = obj.transformLocal2Global(idof);
                dofs = repmat(GlobalDofs',obj.dim.ngaus*obj.dim.nstre,1);
                dofs = dofs(:);
                obj.Btot = obj.Btot + sparse(1:obj.nt,dofs,obj.Bfull(:,idof),obj.nt,obj.ndofGlobal);
            end
        end
        
        
    end
    
end