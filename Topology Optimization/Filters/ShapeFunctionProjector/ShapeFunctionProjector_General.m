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
        
        function xP = project(obj,x)
            if all(x>0)
                xP = zeros(size(x));
            else
                nodalF = ones(size(x));
                obj.unfittedMesh.compute(x);
                s.mesh = obj.unfittedMesh;
                s.type = 'COMPOSITE';
                s = obj.createInteriorParams(s,s.mesh);
                int = Integrator.create(s);
                xP = int.integrateAndSum(nodalF);
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
        
        function cParams = createInteriorParams(obj,cParams,mesh)
            
            if mesh.innerCutMesh.nelem ~= 0
                cParamsInnerCut = obj.createInnerCutParams(mesh);
                cParams.compositeParams{1} = cParamsInnerCut;
                if mesh.innerMesh.nelem ~= 0
                    cParamsInner = obj.createInnerParams(mesh);
                    cParams.compositeParams{end+1} = cParamsInner;
                end
            else
                if mesh.innerMesh.nelem ~= 0
                    cParamsInner = obj.createInnerParams(mesh);
                    cParams.compositeParams{1} = cParamsInner;
                else
                    cParams.compositeParams = cell(0);
                end
            end
        end
        
        function cParams = createInnerParams(obj,mesh)
            cParams.mesh = mesh.innerMesh;
            cParams.type = 'SIMPLE';
            cParams.globalConnec = mesh.globalConnec;
            cParams.npnod = mesh.innerMesh.npnod;
            cParams.backgroundMesh = mesh.meshBackground;
            cParams.innerToBackground = mesh.backgroundFullCells;
        end
        
        function cParams = createInnerCutParams(obj,mesh)
            cParams.mesh = mesh.innerCutMesh;
            cParams.type = 'CutMesh';
            cParams.meshBackground = mesh.meshBackground;
        end
        
    end
    
end