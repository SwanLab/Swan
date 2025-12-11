classdef Data < handle

    properties (Access = public)
        nFeatures
        nSamples
        nLabels        
        Xtrain
        Ytrain       
        Xtest
        Ytest
        Ntest
    end

    properties (Access = private)
        X
        Y
        polynomialOrder
        data
        fileName
        testRatio
        xFeatures
        yFeatures
    end

    methods (Access = public)

        function obj = Data(cParams)            
            obj.init(cParams)
            obj.loadData();
            obj.X = obj.buildModel(obj.X,obj.polynomialOrder);
            obj.splitdata()
            obj.nLabels   = size(obj.Ytrain,2);                        
            obj.nFeatures = size(obj.Xtrain,2);
        end

        function plotdata(self,i,j)
            gscatter(self.Xtrain(:,i),self.Xtrain(:,j),self.Ytrain,'rgbcmyk','*')
            xlabel(['X',num2str(i)]);
            ylabel(['X',num2str(j)]);
        end

        function plotCorrRow(obj,k)
            x = obj.data(:,1:end-1);
            nf = size(x,2);
            for i = 1:nf
                nexttile((i-1)*nf+k)
                if i == k
                    histogram(x(:,i))
                else
                    gscatter(obj.data(:,k),obj.data(:,i),obj.Y,'rgbcmyk','*')
                end
            end
        end

        function plotCorrMatrix(obj)
            x = obj.data(:,1:end-1);
            nf = size(x,2);
            figure            
            t = tiledlayout(nf,nf,'TileSpacing','Compact'); 
            title(t,'Features correlation matrix');
            for i = 1:nf
                obj.plotCorrRow(i);
            end
        end

        function updateHyperparameter(obj,h)
           switch h.type
               case 'testRatio'
                   obj.testRatio = h.value;
                   obj.splitdata()
               case 'polyGrade'
                   obj.polynomialOrder = h.value;
                   obj.buildModel(obj.X,obj.polynomialOrder);
           end
        end


        

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fileName        = cParams.fileName;
            obj.testRatio       = cParams.testRatio;
            obj.polynomialOrder = cParams.polynomialOrder;
            obj.xFeatures       = cParams.xFeatures;
            obj.yFeatures       = cParams.yFeatures;
        end

        function loadData(obj)
            %f = fullfile('../Datasets/',obj.fileName);
            f = fullfile(obj.fileName);

            % Change: use readmatrix to skip header
            obj.data = readmatrix(f);

            % Change: incorporate features to use in cParams vs propmpting
            % user though terminal
            x = obj.data(:, obj.xFeatures);
            y = obj.data(:, obj.yFeatures);

            obj.X = x;
            obj.Y = y;

        end
        

        % function Xful = buildModel(obj)
        %     x  = obj.X;
        %     d  = obj.polynomialOrder;
        %     x1 = x(:,1);
        %     x2 = x(:,2);
        %     cont = 1;
        %     for g = 1:d
        %         for a = 0:g
        %             Xful(:,cont) = x2.^(a).*x1.^(g-a);
        %             cont = cont+1;
        %         end
        %     end
        %     obj.X = Xful;
        % end

  


        function exponents = generateExponents(obj,targetDeg)
            % Initialization of parameters
            exponents = [];
            currentExponents = zeros(1, obj.nFeatures);
            initialIndex = 1;
            
            % Calculation of the possible exponents for the target degree
            exponents = obj.generateExponentsRecursive(targetDeg,initialIndex,currentExponents,exponents);
        end
        
        function exponents = generateExponentsRecursive(obj,targetDeg,currentIndex,currentExponents,exponents)
            % Assignation of exponents for the base case
            if currentIndex == obj.nFeatures
                currentExponents(currentIndex) = targetDeg;
                exponents = [exponents; currentExponents];
            else
                % Recursion to search for the possibile combinations which
                % sum the polynomial degree target
                for i = 0:targetDeg
                    currentExponents(currentIndex) = i;
                    exponents = obj.generateExponentsRecursive(targetDeg - i, currentIndex + 1, currentExponents, exponents);
                end
            end
        end

        function splitdata(obj)
            nD = size(obj.data,1);
            TP = obj.testRatio;
            r = randperm(nD);
            ntest = round(TP/100*nD);
            ntrain = nD - ntest;
            obj.Xtrain = obj.X(r(1:ntrain),:);
            obj.Xtest  = obj.X(r((ntrain + 1):end),:);
            obj.Ytrain = obj.Y(r(1:ntrain),:);
            obj.Ytest  = obj.Y(r((ntrain + 1):end),:);
            obj.Ntest = ntest;
        end
    end


    methods (Access = public, Static)

       function Xful = buildModel(x,d)
            [N,n] = size(x);
            Xful = [];
            for g = 1:d
                E = Data.genExp(n,g);
                for k = 1:size(E,1)
                    c = ones(N,1);
                    for j = 1:n
                        c = c .* x(:,j).^E(k,j);
                    end
                    Xful = [Xful c];
                end
            end
       end

        function E = genExp(n,g)
            if n==1
                E = g;
                return
            end
            E = [];
            for k = 0:g
                S = genExp(n-1,g-k);
                E = [E; k*ones(size(S,1),1) S];
            end
        end       

    end

end