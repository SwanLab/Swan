classdef Preprocess<handle
    % Only reads
    
    %% !! SUSCEPTIBLE OF BEING RENAMED !!
    
    properties
    end
    
    methods(Static)
        function data = readFromGiD(filename)
            
            run(filename)
            data = struct;
            if exist('gidcoord','var')
                coord=gidcoord;
                connec=gidlnods;
            end
            data.xpoints=coord;
            if length(data.xpoints(1,:))==3
                data.xpoints(:,4)=0;
            end
            if any(connec(:,length(connec(1,:))))==0
                data.connectivities=connec(:,1:(length(connec(1,:))-1));
            else
                data.connectivities=connec;
            end
            data.geometry = strjoin(Data_prb(1));
            data.problem_dim = strjoin(Data_prb(3));
            data.problem_type = strjoin(Data_prb(5));
            data.scale = strjoin(Data_prb(6));
            
            if strcmpi(data.problem_type,'elastic')...
               || strcmpi(data.problem_type,'hyperelastic')...
               || strcmpi(data.problem_type,'thermal')
                if exist('dirichlet_data','var')
                    data.dirichlet_data = dirichlet_data;
                elseif exist('lnodes','var')
                    data.dirichlet_data = lnodes;
                else
                end
                
                if exist('pointload_complete','var')
                    data.pointload = pointload_complete;
                else
                    data.pointload = [];
                end
            end
        end
        
        function [fixnodes,forces,boundaryNodes,boundaryElements,Master_slave, sDir, sPL, sPer] = getBC_mechanics(filename)
            run(filename)
            if exist('lnodes','var')
                dirichlet_data=lnodes;
            end
            fixnodes{1,1} = dirichlet_data;
            
            if exist('pointload_complete','var')
                forces = pointload_complete;
            else
                forces = [];
            end
            
            if exist('External_border_nodes','var')
                boundaryNodes= External_border_nodes;
            else
                boundaryNodes = [];
            end

            if exist('External_border_elements','var')
                boundaryElements = External_border_elements;
            else
                boundaryElements = [];
            end            
            
            if ~isempty(boundaryNodes)
                boundaryNodes(:,2)=ones(length(boundaryNodes(:,1)),1);
                boundaryNodes(:,3)=zeros(length(boundaryNodes(:,1)),1);
            end
            
            if ~exist('Master_slave','var')
                Master_slave = [];
            end

            if exist('sDir','var')
                 sDir = sDir;
            else
                 sDir = [];
            end

            if exist('sPL','var')
                 sPL = sPL;
            else
                 sPL = [];
            end

            if exist('sPer','var')
                 sPer = sPer;
            else
                 sPer = [];
            end
            
        end
        
        function [state, velocity, pressure, Vol_force, velocityBC, dtime, finalTime] = getBCFluidsNew(fileName)
            run(fileName)
            if ~exist('dtime', 'var')
                % Steady
                dtime = Inf;
                finalTime = [];
            end
        end

        function [fixnodes,forces,full_dirichlet_data,Master_slave] = getBC_fluids(filename,mesh,geometry,interp)
            run(filename)
            obj = Preprocess;
            full_dirichlet_data = obj.getFullDirichletData(External_border_nodes);

            
            if ~exist('Master_slave','var')
                Master_slave = [];
            end
            
            s.mesh = mesh;
            s.interpolation = interp{1};
            c = ConnecCoordFromInterpAndMesh(s);
            c.compute();
            xpoints = c.coord;
            
            fixnodes_u = obj.getFixedVelocityNodes(xpoints, velocity);
            fixnodes_p = obj.getFixedPressureNodes(pressure);
            forces = obj.getVolumetricForces(Vol_force, mesh, geometry, interp);
            fixnodes{1} = fixnodes_u;
            fixnodes{2} = fixnodes_p;
        end
        
        function forces_adjoint=getBC_adjoint(filename)
            run(filename)
            % for i = 1:numel(sPL)
            %     pl = PointLoad(obj.mesh, sPL{i});
            % end
            forces_adjoint = pointload_adjoint;
        end
    end

    methods (Static, Access = private)
        
        function data = getFullDirichletData(externalBorderNodes)
            fullDirichletData = externalBorderNodes;
            if ~isempty(fullDirichletData)
                fullDirichletData(:,2) = ones(length(fullDirichletData(:,1)),1);
                fullDirichletData(:,3) = zeros(length(fullDirichletData(:,1)),1);
            end
            data = fullDirichletData;
        end

        function nodes = getFixedVelocityNodes(xpoints, velocity)
            nnode   = length(xpoints(:,1));
            
            if (~isempty(velocity))
                ind=1;
                
                for inode = 1: nnode
                    if xpoints(inode,1) ==0  || xpoints(inode,1) == 1 ...
                            || xpoints(inode,2) == 0 || xpoints(inode,2) == 1
                        fixnodes(ind,:) = [inode 1 0];
                        fixnodes(ind+1,:) = [inode 2 0];
                        ind = ind+2;
                        
                    end
                end
                nodes=fixnodes;
                
            end
        end

        function nodes = getFixedPressureNodes(pressure)
            if (~isempty(pressure))
                % if strcmp(geometry(2).interpolation.order,geometry(.order) ==1
                fixnodes_p = pressure;
                %                 else
                %                 ind2=1;
                %
                %                 for i = 1:length(pressure(:,1))
                %
                % %                         [ind,~] = find(interpolation_geometry.T == geometry_fixnodes_p(i,1)); % este para los que tienen minimo un nodo en la pared
                %                    % siguiente linea para los que tienen una cara en la pared
                %                 [ind,~] = find(external_elements(:,2:3) == geometry_fixnodes_p(i,1)); % troba l'element on esta imposada la p a la geometria
                %                 fixnodes_p(ind2: ind2+length(ind)-1,:) = [external_elements(ind,1) ones(length(ind),1)*geometry_fixnodes_p(i,2:3)]; % se li aplica la condicio al primer dels elements
                %                     ind2 = ind2+length(ind);
                %                 end
                %                 fixnodes_p = unique(fixnodes_p,'rows');
                %                 end
                %                 for i = 1:length(obj.fixnodes_p(:,1))
                %                     obj.iD_p(i) = obj.fixnodes_p(i,1)*nunkn_p - nunkn_p + obj.fixnodes_p(i,2);
                %                 end
            end
            nodes = fixnodes_p;
        end

        function forces = getVolumetricForces(Vol_force, mesh, geometry, interp)
            if (~isempty(Vol_force))
                %                 ind=1;
                
                %                 for inode = 1:length(interpolation_variable(1).xpoints(:,1))
                %                     pos_node= num2cell(interpolation_variable(1).xpoints(inode,1:2));
                %                     f = cell2mat(force(pos_node{:}));
                %                     F(ind:ind + length(f)-1,:) = [[inode;inode] [1;2] f];
                %                     ind=ind+length(f);
                %                 end
                geom = geometry(1);
                
                quadrature = Quadrature.set(mesh.type);
                quadrature.computeQuadrature(interp{1}.order);

                geom.computeGeometry(quadrature,interp{1})
                xV = quadrature.posgp;
                xGauss = mesh.computeXgauss(xV);
                nelem = mesh.nelem;
                for ielem = 1:nelem
                    ind=1;
                    for igaus = 1:quadrature.ngaus
                        xG = xGauss(:,igaus,ielem);
                        pos_node= num2cell(xG);
                        f = cell2mat(Vol_force(pos_node{:}));
                        F(:,igaus,ielem) = f;
                        ind=ind+length(f);
                    end
                end
                
                forces = F;
            else
                forces = [];
            end
        end

    end

end