classdef FilterP0 < handle
    
   properties (Access = private)
       levelSet
       levelSet0
       dens0
   end
    
   methods (Access = public)
      
       function obj = FilterP0(ls,d)
           obj.levelSet = ls;
           obj.computeElementalLevelSet(d);
           obj.computeDensityP0();
       end
       
       function d = getDensity(obj)
           d = obj.dens0;
       end
       
   end
   
   methods (Access = private)
       
       function computeElementalLevelSet(obj,d)
            shape = d.shape;
            conec = d.conec;
            quadr = d.quadr;
            ngaus = quadr.ngaus;
            nelem = size(conec,1);
            nnode = size(shape,1);
            
            phiP0 = zeros(nelem,ngaus);
            phi   = obj.levelSet;
            for igaus = 1:ngaus
                for inode = 1:nnode
                    nodes = conec(:,inode);
                    phiN = phi(nodes);
                    phiP0(:,igaus) = phiP0(:,igaus) + shape(inode,igaus)*phiN;
                end
            end
            obj.levelSet0 = phiP0;
       end
       
        function computeDensityP0(obj)
            ls = obj.levelSet0;
            obj.dens0 = obj.computeDensity(ls);
        end
        
   end
   
   methods (Access = private, Static)
        
        function dens = computeDensity(phi)
            dens = 1 - heaviside(phi);
        end
       
   end
   
end