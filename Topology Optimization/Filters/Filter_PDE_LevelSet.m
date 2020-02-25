classdef Filter_PDE_LevelSet < Filter_PDE
    
    properties(Access = private)
        integrator
        unfittedMesh
        domainType
        interp
    end
    
    methods (Access = public)
        
        function obj = Filter_PDE_LevelSet(cParams)
            obj.init(cParams)
            obj.domainType = cParams.domainType;
            
            obj.diffReacProb.preProcess();
            obj.createQuadrature();
            obj.createInterpolation();
            obj.computeGeometry();
            obj.nelem = obj.mesh.nelem;
            obj.npnod = obj.interp.npnod;
            obj.ngaus = obj.quadrature.ngaus;
            obj.Anodal2Gauss = obj.computeA();
            
            %             cParams = SettingsMeshUnfitted(obj.domainType,obj.mesh,obj.interpolation);
            cParams = SettingsMeshUnfitted(obj.domainType,obj.mesh);
            obj.unfittedMesh = UnfittedMesh(cParams);
            % s.mesh = obj.unfittedMesh;
            % s.type = obj.unfittedMesh.unfittedType;
            % obj.integrator = Integrator.create(s);
            obj.disableDelaunayWarning();
        end
        
        function preProcess(obj)
            
        end
        
        function RHS = integrate_L2_function_with_shape_function(obj,x)
            F = ones(size(x));
            RHS = obj.computeRHS(x,F);
        end
        
        function RHS = integrate_function_along_facets(obj,x,F)
            RHS = obj.computeRHS(x,F);
        end
        
        function fInt = computeRHS(obj,ls,fNodes)
            if all(ls>0)
                fInt = zeros(size(ls));
            else
                obj.unfittedMesh.compute(ls);                
                s.mesh = obj.unfittedMesh;
                s.type = 'COMPOSITE';
                s = obj.createInteriorParams(s,s.mesh);
                int = Integrator.create(s);
                fInt = int.integrateAndSum(fNodes);
            end
        end
        
    end
    
    methods (Access = private)
        
        function createQuadrature(obj)
            obj.quadrature = Quadrature.set(obj.mesh.geometryType);
            obj.quadrature.computeQuadrature(obj.quadratureOrder);
        end
        
        function createInterpolation(obj)
            obj.interp = Interpolation.create(obj.mesh,'LINEAR');
        end
        
        function computeGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry = Geometry.create(s);
            obj.geometry.computeGeometry(obj.quadrature,obj.interp);
        end
        
        
        function disableDelaunayWarning(obj)
            MSGID = 'MATLAB:delaunayTriangulation:DupPtsWarnId';
            warning('off', MSGID)
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