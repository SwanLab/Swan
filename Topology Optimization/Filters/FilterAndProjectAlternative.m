classdef FilterAndProjectAlternative < Filter
    properties (Access = public)
        filteredField
        projectedField
    end
    properties (Access = private)
        eta
        beta
        filter
        projector
        weight
        weightSum
        
    end
    
    methods (Access = public)
        function obj = FilterAndProjectAlternative(cParams)
            obj.init(cParams);
            obj.createWeight()
            obj.defineProjectorSettings(cParams)
            obj.createFilter();
        end

        function x_reg_Elements = getP0fromP1(obj,x)
            obj.computeFilter(x);
            obj.createProjector(); 
            obj.computeProjector();
            x_reg = obj.projector.projectedField;
            x_reg_Elements = FilterAndProjectAlternative.nodesToElements(x_reg,obj.mesh);
        end

        function x_reg = getP1fromP0(obj,dC_dRhoBar)
            s.beta = obj.beta;
            s.eta  = obj.eta;
            s.filteredField = obj.filteredField;
            B = ProjectedFieldFilteredFieldDerivator(s);
            B.compute();
            dRhoBar_dRhoHat = B.derivatedProjectedField;

            s.H = obj.weight;
            s.Hs = obj.weightSum;
            B = FilteredFieldFieldDerivator(s);
            B.compute();
            dRhoHat_dRho = B.derivedFilteredField;
            x_reg = FilterAndProjectAlternative.elementsToNodes(dC_dRhoBar,obj.mesh).*(dRhoHat_dRho*dRhoBar_dRhoHat);
        end

    end
    methods (Access = private)
        function createFilter(obj)
            s.filterParameters.H = obj.weight;
            s.filterParameters.Hs = obj.weightSum;
            obj.filter  = FilterComputer(s);
        end
        function defineProjectorSettings(obj,cParams)
            obj.eta  = cParams.femSettings.eta;
            obj.beta = cParams.femSettings.beta;
        end 
        function computeFilter(obj,density)
            obj.filteredField = obj.filter.compute(density);
        end
        function createProjector(obj)
            s.beta = obj.beta;
            s.eta = obj.eta;
            s.filteredField = obj.filteredField ;
            obj.projector = FieldProjector(s);
        end
        function computeProjector(obj)
            obj.projector.compute();
        end 
        function createWeight(obj)
            [XnodesNumber,YnodesNumber] = FilterAndProjectAlternative.countNodesXY(obj.mesh.coord);
            s.Xnumber =XnodesNumber;
            s.Ynumber = YnodesNumber;
            s.minimunInfluenceRadios = 1;
            B = weightFilterComputer(s);
            B.compute;
            obj.weight = B.H;
            obj.weightSum = B.Hs;
        end  
    end
    methods(Static, Access = protected)
        function [XnodesNumber,YnodesNumber] = countNodesXY(coord)
            XnodesNumber = numel(unique(coord(:,1)));
            YnodesNumber = numel(unique(coord(:,2)));
        end
        function  valueInElements = nodesToElements(valueInNodes,mesh)
            valueInElements = zeros(size(mesh.connec));
            connectivityMat = mesh.connec;
            for i = 1:mesh.nelem
                valueInElements(i,:) = valueInNodes(connectivityMat(i,:))';
            end
        end
        function valueInNodes = elementsToNodes(valueInElements,mesh)
            valueInNodes = zeros(mesh.nnodes,1);
            for i = 1:mesh.nelem
                nodes = mesh.connec(i, :);
                valueInNodes(nodes) = valueInNodes(nodes) + valueInElements(i, :)';
            end
        end
    end
end