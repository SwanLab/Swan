classdef ShapeFunctionProjector_General < ShapeFunctionProjector
    
    properties (Access = private)
        unfittedMesh
        domainType
        interpolation
        quadrature
    end
    
    methods (Access = public)
        
        function obj = ShapeFunctionProjector_General(cParams)
            obj.init(cParams)
            obj.domainType = cParams.domainType;
            obj.quadrature = cParams.quadrature;
            obj.createInterpolation();
            obj.createUnfittedMesh();
        end
        
        function fInt = project(obj,ls)
            if all(ls>0)
                fInt = zeros(size(ls));
            else
                fNodes = ones(size(ls));
                obj.unfittedMesh.compute(ls);
%                 close all
%                 figure(1)
%                 obj.unfittedMesh.plot();
%                 drawnow
%                 pause(1)
                fInt = obj.unfittedMesh.integrateNodalFunction(fNodes);
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function createInterpolation(obj)
            obj.interpolation = Interpolation.create(obj.mesh,'LINEAR');
            obj.interpolation.computeShapeDeriv(obj.quadrature.posgp)
        end
        
        function createUnfittedMesh(obj)
            s.unfittedType = obj.domainType;
            s.meshBackground = obj.mesh;
            s.interpolationBackground = obj.interpolation;
            cParams = SettingsMeshUnfitted(s);
            obj.unfittedMesh = UnfittedMesh(cParams);
        end
        
    end
    
end