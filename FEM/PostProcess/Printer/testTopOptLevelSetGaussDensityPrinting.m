classdef testTopOptLevelSetGaussDensityPrinting < testTopOptPrinting
    
    properties (Access = protected)
        testName = 'test_micro2';  
        fileOutputName = 'testTopOptLevelSetGaussDensityPrinting';
        postProcessor = 'LevelSetGaussDensity';
    end
 
    properties (Access = private)
        ls
        lsP0 
        densP0
    end
    
    methods (Access = protected)
       
        function computeFields(obj)
            obj.computePhiP0();
            obj.computeDensityP0();
            obj.fields.density  = obj.densP0;
            obj.fields.levelSet = obj.ls;
        end
        
        function computePhiP0(obj)            
            conec = obj.topOpt.mesh.connec;
            phyPr = obj.topOpt.cost.ShapeFuncs{1}.physicalProblem;            
            shape = phyPr.element.interpolation_u.shape;
            quadr = phyPr.element.quadrature;
            ngaus = quadr.ngaus;
            nelem = size(conec,1);
            nnode = size(shape,1);

            phiP0 = zeros(ngaus,nelem);
            phi   = obj.topOpt.x;
            
            for igaus = 1:ngaus
                for inode = 1:nnode 
                    nodes = conec(:,inode);
                    phiN(1,:) = phi(nodes);
                    phiP0(igaus,:) = phiP0(igaus,:) + shape(inode,igaus)*phiN;
                end
            end            
            obj.lsP0 = phiP0;
            obj.ls   = phi;
        end
        
        function computeDensityP0(obj)
            obj.densP0 = 1 - heaviside(obj.lsP0);            
        end
        
        function d = createPostProcessDataBaseStructre(obj)
            fem        = obj.topOpt.cost.ShapeFuncs{1}.physicalProblem;
            dI.mesh    = obj.topOpt.mesh;
            dI.fields  = obj.fields;
            dI.outName = obj.fileOutputName;
            dI.quad    = fem.element.quadrature;
            dI.iter    = obj.iter;
            hasGaussData = true;
            ps = PostProcessDataBaseCreator.create(hasGaussData,dI);
            d = ps.getValue();              
        end
  
        
    end
end
