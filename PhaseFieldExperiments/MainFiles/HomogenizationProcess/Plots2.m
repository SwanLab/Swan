matType{1} = load('CircleMicroDamageArea.mat');
matType{2} = load('CircleMicroDamagePerimeter.mat');
matType{3} = load('SquareMicroDamageArea.mat');
matType{4} = load('SquareMicroDamagePerimeter.mat');
matType{5} = load('IsoMicroDamage.mat');
matType{6} = load('CrackMicroDamageCubic.mat');
matType{6}.mat = matType{6}.mat*210;

% Change of variable
%matType{6}.phi  = matType{6}.holeParam{1};
%matType{6}.phi = matType{6}.holeParam{1}.^2;
matType{6}.phi = matType{6}.holeParam{1}.^3;

%% Polynomial fraction fitting
close all
x = matType{6}.phi;
y = squeeze(matType{6}.mat(3,3,:));
p0 = [1 1 1 1 1 1 1 1];
yp = @(p) (p(1).*x.^3 + p(2).*x.^2 + p(3).*x + p(4))/(p(5).*x.^3 + p(6).*x.^2 + p(7).*x + p(8));

objective = @(p) sum(((yp(p)-y)./y).^2);
popt = fmincon(objective,p0);

disp("Initial objective: " + num2str(objective(p0)));
disp("Final objective: " + num2str(objective(p0)));
plot(x,y,'ro')
hold on
plot(x,yp(popt),'gs')
legend('measured','optimal')

%% Include pints
% %% Include final points
% matType{6}.mat(:,:,end+1) = [matType{6}.mat(1,1,end), 0 , 0;
%                                         0           , 0 , 0;
%                                         0           , 0 , 0];
% matType{6}.phi(end+1) = 1;
% matType{6}.mat(1,1,end) = matType{6}.mat(1,1,end-1);
%
% %% Include more final points %%
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
     [funMat(:,:,j),dfunMat(:,:,j),ddfunMat(:,:,j)] = computeFittingHomogenization(matType{j},isMat6,4);
 end

 %% Run Plots
 Finalplots;




