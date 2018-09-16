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
    end
    
    methods
        
        function obj = KGeneratorWithfullStoredB(Bfull,dim,connectivities,Cmat,dvolum)
            
            obj.dofsPerElement   = dim.ndofPerElement;
            obj.nodesInElement   = reshape(repmat(1:dim.nnode,dim.nunkn,1),1,[]);
            obj.VectorDimensions = repmat(1:dim.nunkn,1,dim.nnode);
            obj.connectivities   = connectivities;
            
            obj.ndofGlobal = max(max(connectivities))*dim.nunkn;
            nt = dim.ngaus*dim.nelem*dim.nstre;
            obj.Btot = sparse(nt,obj.ndofGlobal);
            obj.Bfull = Bfull;
            obj.dim = dim; 
            obj.nunkn = dim.nunkn;
            obj.computeBtot();
            obj.computeCmatBlockDiagonal(Cmat,dvolum);
            obj.computeStiffnes()
        end
        
        function  computeStiffnes(obj)

            B = obj.Btot;
            CB = obj.CmatTot*B;
            obj.K = B'*CB;            
        end
        
        
        function computeCmatBlockDiagonal(obj,Cmat,dvolum)
            nt = obj.dim.ngaus*obj.dim.nelem*obj.dim.nstre;
            obj.CmatTot = sparse(nt,nt);
            for istre = 1:obj.dim.nstre
                for jstre = 1:obj.dim.nstre
                    for igaus = 1:obj.dim.ngaus
                        posI = (istre)+(obj.dim.nstre)*(igaus-1) : obj.dim.ngaus*obj.dim.nstre : nt ;
                        posJ = (jstre)+(obj.dim.nstre)*(igaus-1) : obj.dim.ngaus*obj.dim.nstre : nt ;
                        
                        Ct = squeeze(Cmat(istre,jstre,:)).*dvolum(:,igaus);                        
                        obj.CmatTot = obj.CmatTot + sparse(posI,posJ,Ct,nt,nt);
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
            
            for idof=1:obj.dofsPerElement
                GlobalDofs = transformLocal2Global(obj,idof);
                
                dofs = repmat(GlobalDofs',obj.dim.ngaus*obj.dim.nstre,1);
                dofs = dofs(:);
                nt = obj.dim.ngaus*obj.dim.nelem*obj.dim.nstre;
                obj.Btot = obj.Btot + sparse(1:nt,dofs,obj.Bfull(:,idof),nt,obj.ndofGlobal);
            end
        end
        
        
    end
    
end