classdef FilterP0 < handle
    
   properties (Access = private)
       levelSet
       levelSet0
       phyProb       
       dens0       
   end
    
   methods (Access = public)
      
       function obj = FilterP0(ls,phy)
           obj.levelSet = ls;
           obj.phyProb = phy;           
           obj.computeElementalLevelSet();
           obj.computeDensityP0();
       end
       
       function d = getDens0(obj)
           d = obj.dens0;
       end
       
   end
   
   methods (Access = private)
       
       function computeElementalLevelSet(obj)
            phyPr = obj.phyProb();
            shape = phyPr.element.interpolation_u.shape;
            conec = phyPr.geometry.interpolation.T;            
            quadr = phyPr.element.quadrature;
            ngaus = quadr.ngaus;
            nelem = size(conec,1);
            nnode = size(shape,1);
            
            phiP0 = zeros(ngaus,nelem);
            phi   = obj.levelSet;            
            for igaus = 1:ngaus
                for inode = 1:nnode
                    nodes = conec(:,inode);
                    phiN(1,:) = phi(nodes);
                    phiP0(igaus,:) = phiP0(igaus,:) + shape(inode,igaus)*phiN;
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