classdef NumericalRectangleTypeHomogenizer < NumericalHomogenizer
    
    properties (Access = protected)
        nodalLevelSet        
    end
    
    properties (Access = private)
        m1
        m2
    end
    
      methods (Access = protected)
        
        function compute(obj,fileName,print,m1,m2,iter)
            obj.m1 = m1;
            obj.m2 = m2;
            obj.init(fileName,print,iter);
            obj.generateMicroProblem();
            obj.computeHomogenizedVariables();
            obj.createDensityPrinter()
            obj.print()
        end
        
    end
    
    
    methods (Access = protected)
        
        function createDensity(obj)
            obj.createLevelSet(obj.m1,obj.m2)                                    
            shape   = obj.microProblem.element.interpolation_u.shape; 
            conec   = obj.microProblem.geometry.interpolation.T;

            nelem = size(conec,1);
            nnode = size(shape,1);
            quadr = obj.microProblem.element.quadrature;
            ngaus = quadr.ngaus;            

            lsGaus = zeros(ngaus,nelem);
            dens = zeros(ngaus,nelem);
            lsNodes = obj.nodalLevelSet;
            for igaus = 1:ngaus
                for inode =  1:nnode
                    nodes = conec(:,inode);
                    ls  = lsNodes(nodes);
                    lsGaus(igaus,:) = lsGaus(igaus,:) + shape(inode,igaus)*ls';
                end
                dens(igaus,lsGaus(igaus,:)<0) = 1;
            end            
            obj.density = dens';
            
        end
        
    end
    
    methods (Access = protected, Abstract)
       createLevelSet(obj,m1,m2) 
    end
    
end

