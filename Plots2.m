matType{1} = load('CircleMicroDamageArea.mat');
matType{2} = load('CircleMicroDamagePerimeter.mat');
matType{3} = load('SquareMicroDamageArea.mat');
matType{4} = load('SquareMicroDamagePerimeter.mat');
matType{5} = load('IsoMicroDamage.mat');
matType{6} = load('CrackMicroDamage100.mat');
matType{6}.mat = matType{6}.mat*210;

%% Change of variable
matType{6}.phi  = matType{6}.holeParam{1};
%matType{6}.phi = matType{6}.holeParam{1}.^2;
%matType{6}.phi = matType{6}.holeParam{1}.^3;

%% Include last point
% matType{6}.mat(:,:,end+1) = [matType{6}.mat(1,1,end), 0 , 0;
%                                         0           , 0 , 0;
%                                         0           , 0 , 0];
% matType{6}.phi(end+1) = 1;
% matType{6}.mat(1,1,end) = matType{6}.mat(1,1,end-1);
% 
% %% Include more final points
% lastSteps = 2;
% nPoints = 20;
% lastPointsMat = matType{6}.mat(:,:,(end-lastSteps):end);
% lastPointsPhi = matType{6}.phi((end-lastSteps):end);
% mat = zeros(3,3,lastSteps*(nPoints-1)+1);
% phi = zeros(1,lastSteps*(nPoints-1)+1);
% for k=1:lastSteps
%     if k==1 
%         extra = 1; 
%     else 
%         extra=0;
%     end
% 
%     for i=1:3
%         for j=1:3
%             start  = lastPointsMat(i,j,k);
%             finish = lastPointsMat(i,j,k+1);
%             res = linspace(start,finish,nPoints);
%             mat(i,j,((k-1)*nPoints+extra:k*nPoints-1)) = res(1:end-extra);
%         end
%     end
%     start  = lastPointsPhi(k);
%     finish = lastPointsPhi(k+1);
%     res = linspace(start,finish,nPoints);
%     phi(((k-1)*nPoints+extra:k*nPoints-1)) = res(1:end-extra);
% end
% mat(:,:,end) = lastPointsMat(:,:,end);
% phi(end) = lastPointsPhi(end);
% 
% matType{6}.mat = cat(3,matType{6}.mat(:,:,1:(end-lastSteps-1)),mat);
% matType{6}.phi = [matType{6}.phi(1:(end-lastSteps-1)),phi];

%% Fitting
 for j=1:length(matType)
     isMat6 = false;
     if j==6 isMat6 = true; end
     [funMat(:,:,j),dfunMat(:,:,j),ddfunMat(:,:,j)] = computeFunctionsAndDerivatives(matType{j},isMat6);
 end

 %% Run Plots
 Finalplots;


 %% Functions 
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
                %coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0],[squeeze(y(i,j,1))]);
                %coeffs = polyfix(x,squeeze(y(i,j,:)),9,[1],[squeeze(y(i,j,end))]);
                coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0,1],[squeeze(y(i,j,1)),squeeze(y(i,j,end))]);
            else
                %coeffs = polyfit(x,squeeze(y(i,j,:)),9);
                %coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0],[squeeze(y(i,j,1))]);
                %coeffs = polyfix(x,squeeze(y(i,j,:)),9,[1],[0]);
                coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0,1],[squeeze(y(i,j,1)),0]);
            end
            fun{i,j} = poly2sym(coeffs);
            dfun{i,j} = diff(fun{i,j});
            ddfun{i,j} = diff(dfun{i,j});
        end
    end
 end

