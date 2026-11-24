% % classdef LatticeDehomogenizer < handle
% % 
% %     properties (Access = private)
% %         design
% %         epsilon
% %         elementsPerCell
% %         rhoVoid
% %         mesh
% %         fields
% %         levelSet
% %         uMesh
% %     end
% % 
% %     methods (Access = public)
% % 
% %         function obj = LatticeDehomogenizer(cParams)
% %             obj.init(cParams);
% %             obj.createMesh();
% %             obj.interpolateFields();
% %             obj.createLevelSet();
% %             obj.createUnfittedMesh();
% %         end
% % 
% %         function plot(obj)
% %             figure
% %             obj.uMesh.plot();
% %             axis equal off
% %         end
% % 
% %         function v = computeVolumeFraction(obj)
% %             chi      = CharacteristicFunction.create(obj.uMesh);
% %             sF.trial = LagrangianFunction.create(obj.mesh,1,'P1');
% %             sF.mesh  = obj.mesh;
% %             dens     = FilterLump(sF).compute(chi,2);
% %             v        = Integrator.compute(dens,obj.mesh,2)/obj.mesh.computeVolume();
% %         end
% % 
% %         function v = computeTargetVolumeFraction(obj)
% %             rho = obj.fields.rho;
% %             rho(rho < obj.rhoVoid) = 0;
% %             f = LagrangianFunction.create(obj.mesh,1,'P1');
% %             f.setFValues(rho);
% %             v = Integrator.compute(f,obj.mesh,2)/obj.mesh.computeVolume();
% %         end
% % 
% %         function ls = getLevelSet(obj)
% %             ls = obj.levelSet;
% %         end
% % 
% %         function m = getMesh(obj)
% %             m = obj.mesh;
% %         end
% % 
% %     end
% % 
% %     methods (Access = private)
% % 
% %         function init(obj,cParams)
% %             obj.design          = load(cParams.designFile);
% %             obj.epsilon         = cParams.epsilon;
% %             obj.elementsPerCell = cParams.elementsPerCell;
% %             obj.rhoVoid         = 1 - obj.design.mU^2;
% %         end
% % 
% %         function createMesh(obj)
% %             coord      = obj.design.coord;
% %             domainSize = max(coord) - min(coord);
% %             nCells     = round(domainSize/obj.epsilon);
% %             nDiv       = nCells*obj.elementsPerCell;
% %             obj.mesh   = TriangleMesh(domainSize(1),domainSize(2),nDiv(1),nDiv(2));
% %         end
% % 
% %         function interpolateFields(obj)
% %             coord = obj.design.coord;
% %             x     = obj.mesh.coord;
% %             names = {'rho','a','xi'};
% %             for k = 1:numel(names)
% %                 F = scatteredInterpolant(coord(:,1),coord(:,2),obj.design.(names{k}),'linear','nearest');
% %                 obj.fields.(names{k}) = F(x(:,1),x(:,2));
% %             end
% %         end
% % 
% %         function createLevelSet(obj)
% %             x   = obj.mesh.coord;
% %             rho = obj.fields.rho;
% %             a   = obj.fields.a;
% %             xi  = obj.fields.xi;
% %             cell = x/obj.epsilon;
% %             cell = cell - floor(cell) - 0.5;
% %             k    = 1./sqrt(1-a.^2);
% %             phi1 = k.*(cell(:,1) - a.*cell(:,2));
% %             phi2 = k.*(cell(:,2) - a.*cell(:,1));
% %             m0   = sqrt(1-rho);
% %             r    = sqrt(tan(xi));
% %             m1   = m0.*r;
% %             m2   = m0./r;
% %             psi  = max(abs(phi1)./(m1/2),abs(phi2)./(m2/2)) - 1;
% %             psi(rho < obj.rhoVoid) = -1;
% %             obj.levelSet = LagrangianFunction.create(obj.mesh,1,'P1');
% %             obj.levelSet.setFValues(psi);
% %         end
% % 
% %         function createUnfittedMesh(obj)
% %             s.backgroundMesh = obj.mesh;
% %             s.boundaryMesh   = obj.mesh.createBoundaryMesh();
% %             obj.uMesh = UnfittedMesh(s);
% %             obj.uMesh.compute(-obj.levelSet.fValues);
% %         end
% % 
% %     end
% % 
% % end
% classdef LatticeDehomogenizer < handle
% 
%     properties (Access = private)
%         design
%         epsilon
%         elementsPerCell
%         rhoVoid
%         mesh
%         fields
%         levelSet
%         uMesh
%     end
% 
%     methods (Access = public)
% 
%         function obj = LatticeDehomogenizer(cParams)
%             obj.init(cParams);
%             obj.createMesh();
%             obj.interpolateFields();
%             obj.createLevelSet();
%             obj.createUnfittedMesh();
%         end
% 
%         function plot(obj)
%             figure
%             obj.uMesh.plot();
%             axis equal off
%         end
% 
%         function v = computeVolumeFraction(obj)
%             chi      = CharacteristicFunction.create(obj.uMesh);
%             sF.trial = LagrangianFunction.create(obj.mesh,1,'P1');
%             sF.mesh  = obj.mesh;
%             dens     = FilterLump(sF).compute(chi,2);
%             v        = Integrator.compute(dens,obj.mesh,2)/obj.mesh.computeVolume();
%         end
% 
%         function v = computeTargetVolumeFraction(obj)
%             rho = obj.fields.rho;
%             rho(rho < obj.rhoVoid) = 0;
%             f = LagrangianFunction.create(obj.mesh,1,'P1');
%             f.setFValues(rho);
%             v = Integrator.compute(f,obj.mesh,2)/obj.mesh.computeVolume();
%         end
% 
%         function ls = getLevelSet(obj)
%             ls = obj.levelSet;
%         end
% 
%         function m = getMesh(obj)
%             m = obj.mesh;
%         end
% 
%     end
% 
%     methods (Access = private)
% 
%         function init(obj,cParams)
%             obj.design          = load(cParams.designFile);
%             obj.epsilon         = cParams.epsilon;
%             obj.elementsPerCell = cParams.elementsPerCell;
%             obj.rhoVoid         = 1 - obj.design.mU^2;
%         end
% 
%         function createMesh(obj)
%             coord      = obj.design.coord;
%             domainSize = max(coord) - min(coord);
%             nCells     = round(domainSize/obj.epsilon);
%             nDiv       = nCells*obj.elementsPerCell;
%             obj.mesh   = TriangleMesh(domainSize(1),domainSize(2),nDiv(1),nDiv(2));
%         end
% 
%         function interpolateFields(obj)
%             coord  = obj.design.coord;
%             connec = obj.design.connec;
%             x      = obj.mesh.coord;
%             TR     = triangulation(connec,coord);
%             [element,baricentric] = pointLocation(TR,x);
%             inside  = ~isnan(element);
%             nearest = nearestNeighbor(TR,x(~inside,:));
%             names   = {'rho','a','xi'};
%             for k = 1:numel(names)
%                 nodal = obj.design.(names{k});
%                 value = zeros(size(x,1),1);
%                 value(inside)  = sum(baricentric(inside,:).*nodal(connec(element(inside),:)),2);
%                 value(~inside) = nodal(nearest);
%                 obj.fields.(names{k}) = value;
%             end
%         end
% 
%         function createLevelSet(obj)
%             x   = obj.mesh.coord;
%             rho = obj.fields.rho;
%             a   = obj.fields.a;
%             xi  = obj.fields.xi;
%             cell = x/obj.epsilon;
%             cell = cell - floor(cell) - 0.5;
%             k    = 1./sqrt(1-a.^2);
%             phi1 = k.*(cell(:,1) - a.*cell(:,2));
%             phi2 = k.*(cell(:,2) - a.*cell(:,1));
%             m0   = sqrt(1-rho);
%             r    = sqrt(tan(xi));
%             m1   = m0.*r;
%             m2   = m0./r;
%             psi  = max(abs(phi1)./(m1/2),abs(phi2)./(m2/2)) - 1;
%             psi(rho < obj.rhoVoid) = -1;
%             obj.levelSet = LagrangianFunction.create(obj.mesh,1,'P1');
%             obj.levelSet.setFValues(psi);
%         end
% 
%         function createUnfittedMesh(obj)
%             s.backgroundMesh = obj.mesh;
%             s.boundaryMesh   = obj.mesh.createBoundaryMesh();
%             obj.uMesh = UnfittedMesh(s);
%             obj.uMesh.compute(-obj.levelSet.fValues);
%         end
% 
%     end
% 
% end

classdef LatticeDehomogenizer < handle
    properties (Access = private)
        design
        epsilon
        elementsPerCell
        minimumElementsPerWall
        rhoVoid
        mesh
        fields
        sides
        transformation
        parameterNames
        levelSet
        uMesh
    end
    methods (Access = public)
        function obj = LatticeDehomogenizer(cParams)
            obj.init(cParams);
            obj.createMesh();
            obj.interpolateFields();
            obj.createLevelSet();
            obj.createUnfittedMesh();
        end
        function plot(obj)
            figure
            obj.uMesh.plot();
            axis equal off
        end
        function v = computeVolumeFraction(obj)
            chi      = CharacteristicFunction.create(obj.uMesh);
            sF.trial = LagrangianFunction.create(obj.mesh,1,'P1');
            sF.mesh  = obj.mesh;
            dens     = FilterLump(sF).compute(chi,2);
            v        = Integrator.compute(dens,obj.mesh,2)/obj.mesh.computeVolume();
        end
        function v = computeTargetVolumeFraction(obj)
            rho = obj.fields.rho;
            rho(rho < obj.rhoVoid) = 0;
            f = LagrangianFunction.create(obj.mesh,1,'P1');
            f.setFValues(rho);
            v = Integrator.compute(f,obj.mesh,2)/obj.mesh.computeVolume();
        end
        function v = computeRepresentableVolumeFraction(obj)
            rho = 1 - obj.sides(:,1).*obj.sides(:,2);
            f = LagrangianFunction.create(obj.mesh,1,'P1');
            f.setFValues(rho);
            v = Integrator.compute(f,obj.mesh,2)/obj.mesh.computeVolume();
        end
        function ls = getLevelSet(obj)
            ls = obj.levelSet;
        end
        function m = getMesh(obj)
            m = obj.mesh;
        end
    end
    methods (Access = private)
        function init(obj,cParams)
            obj.design          = load(cParams.designFile);
            obj.epsilon         = cParams.epsilon;
            obj.elementsPerCell = cParams.elementsPerCell;
            obj.rhoVoid         = 1 - obj.design.mU^2;
            obj.minimumElementsPerWall = 1;
            if isfield(cParams,'minimumElementsPerWall')
                obj.minimumElementsPerWall = cParams.minimumElementsPerWall;
            end
            obj.transformation = BoostTransformation();
            if isfield(cParams,'transformation')
                obj.transformation = cParams.transformation;
            end
            obj.parameterNames = {'a'};
            if isfield(cParams,'parameterNames')
                obj.parameterNames = cParams.parameterNames;
            end
        end
        function createMesh(obj)
            coord      = obj.design.coord;
            domainSize = max(coord) - min(coord);
            nCells     = round(domainSize/obj.epsilon);
            nDiv       = nCells*obj.elementsPerCell;
            obj.mesh   = TriangleMesh(domainSize(1),domainSize(2),nDiv(1),nDiv(2));
        end
        function interpolateFields(obj)
            coord  = obj.design.coord;
            connec = obj.design.connec;
            x      = obj.mesh.coord;
            TR     = triangulation(connec,coord);
            [element,baricentric] = pointLocation(TR,x);
            inside  = ~isnan(element);
            nearest = nearestNeighbor(TR,x(~inside,:));
            % names   = {'rho','a','xi'};
            names   = [{'rho','xi'},obj.parameterNames];
            for k = 1:numel(names)
                nodal = obj.design.(names{k});
                value = zeros(size(x,1),1);
                value(inside)  = sum(baricentric(inside,:).*nodal(connec(element(inside),:)),2);
                value(~inside) = nodal(nearest);
                obj.fields.(names{k}) = value;
            end
        end
        function createLevelSet(obj)
            x   = obj.mesh.coord;
            rho = obj.fields.rho;
            % a   = obj.fields.a;
            nParams = numel(obj.parameterNames);
            params  = zeros(size(x,1),nParams);
            for p = 1:nParams
                params(:,p) = obj.fields.(obj.parameterNames{p});
            end
            xi  = obj.fields.xi;
            cell = x/obj.epsilon;
            cell = cell - floor(cell) - 0.5;
            % k    = 1./sqrt(1-a.^2);
            % phi1 = k.*(cell(:,1) - a.*cell(:,2));
            % phi2 = k.*(cell(:,2) - a.*cell(:,1));
            [values,~,index] = unique(params,'rows');
            Tinv = zeros(size(values,1),2,2);
            for i = 1:size(values,1)
                Tinv(i,:,:) = inv(obj.transformation.compute(values(i,:)));
            end
            phi1 = Tinv(index,1,1).*cell(:,1) + Tinv(index,1,2).*cell(:,2);
            phi2 = Tinv(index,2,1).*cell(:,1) + Tinv(index,2,2).*cell(:,2);
            m0   = sqrt(1-rho);
            r    = sqrt(tan(xi));
            m1   = m0.*r;
            m2   = m0./r;
            % psi  = max(abs(phi1)./(m1/2),abs(phi2)./(m2/2)) - 1;
            % psi(rho < obj.rhoVoid) = -1;
            isWall1 = obj.isWallResolvable(m1);
            isWall2 = obj.isWallResolvable(m2);
            term1   = zeros(size(m1));
            term2   = zeros(size(m2));
            term1(isWall1) = abs(phi1(isWall1))./(m1(isWall1)/2);
            term2(isWall2) = abs(phi2(isWall2))./(m2(isWall2)/2);
            psi = max(term1,term2) - 1;
            m1(~isWall1) = 1;
            m2(~isWall2) = 1;
            obj.sides    = [m1,m2];
            obj.levelSet = LagrangianFunction.create(obj.mesh,1,'P1');
            obj.levelSet.setFValues(psi);
        end
        function is = isWallResolvable(obj,m)
            is = (1 - m)*obj.elementsPerCell >= obj.minimumElementsPerWall;
        end
        function createUnfittedMesh(obj)
            s.backgroundMesh = obj.mesh;
            s.boundaryMesh   = obj.mesh.createBoundaryMesh();
            obj.uMesh = UnfittedMesh(s);
            obj.uMesh.compute(-obj.levelSet.fValues);
        end
    end
end