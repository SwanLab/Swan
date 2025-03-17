matType{1} = load('CircleMicroDamageArea.mat');
matType{2} = load('CircleMicroDamagePerimeter.mat');
matType{3} = load('SquareMicroDamageArea.mat');
matType{4} = load('SquareMicroDamagePerimeter.mat');
matType{5} = load('IsoMicroDamage.mat');
matType{6} = load('HorizontalCrackMicroDamageArea.mat');
matType{6}.mat = matType{6}.mat*210;

mat11 = squeeze(matType{6}.mat(1,1,:))';
mat12 = squeeze(matType{6}.mat(1,2,:))';
phi   = squeeze(matType{6}.phi);
curveFitter(phi,mat12);

%% Change of variable
matType{6}.phi = matType{6}.holeParam{1};
%matType{6}.phi = matType{6}.holeParam{1}.^2;
%matType{6}.phi = matType{6}.holeParam{1}.^3;

%% Include last point
% matType{6}.mat(:,:,end+1) = [matType{6}.mat(1,1,end), 0 , 0;
%                                         0           , 0 , 0;
%                                         0           , 0 , 0];
% matType{6}.phi(end+1) = 1;

%% Include final points
% lastMat = matType{6}.mat(:,:,end);
% lastPhi = matType{6}.phi(end);
% nPoints = 20;
% points = zeros(3,3,nPoints);
% for i=1:3
%     for j=1:3
%         if i==1 && j==1
%             points(i,j,:) = lastMat(i,j);
%         else
%             points(i,j,:) = linspace(lastMat(i,j),0,nPoints); 
%         end
%     end
% end
% phiPoints = linspace(lastPhi,1,nPoints);
% 
% matType{6}.mat = cat(3,matType{6}.mat,points);
% matType{6}.phi = [matType{6}.phi,phiPoints];

%% Run Plots
Finalplots;

%% Fitting
 for i=1:length(matType)
     isMat6 = false;
     if i==6 isMat6 = true; end
     [funMat(:,:,i),dfunMat(:,:,i),ddfunMat(:,:,i)] = computeFunctionsAndDerivatives(matType{i},isMat6);
 end

 function [fun,dfun,ddfun] = computeFunctionsAndDerivatives(cParams,isMat6)
    x = cParams.phi;
    x = reshape(x,length(x),[]);
    y = cParams.mat;
    
    fun   = cell(3,3);
    dfun  = cell(3,3);
    ddfun = cell(3,3);
    for i=1:3
        for j=1:3
            if i==1 && j==1 && isMat6
                %coeffs = polyfit(x,squeeze(y(i,j,:)),9);
                coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0],[squeeze(y(i,j,1))]);
                %coeffs = polyfix(x,squeeze(y(i,j,:)),9,[1],[squeeze(y(i,j,end))]);
                %coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0,1],[squeeze(y(i,j,1)),squeeze(y(i,j,end))]);
            else
                %coeffs = polyfit(x,squeeze(y(i,j,:)),9);
                coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0],[squeeze(y(i,j,1))]);
                %coeffs = polyfix(x,squeeze(y(i,j,:)),9,[1],[0]);
                %coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0,1],[squeeze(y(i,j,1)),0]);
            end
            fun{i,j} = poly2sym(coeffs);
            dfun{i,j} = diff(fun{i,j});
            ddfun{i,j} = diff(dfun{i,j});
        end
    end
 end