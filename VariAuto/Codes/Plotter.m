classdef Plotter < handle

    properties (Access = private)
        data
        neuronsPerLayer
        network
    end

    methods (Access = public)

        function self = Plotter(init)
            self.data = init.data;
            self.network = init;
            self.neuronsPerLayer = init.network.neuronsPerLayer;
        end

        function plotBoundary(self,type) 
           X = self.data.Xtrain;
           nF = size(X,2);
           nPL = self.neuronsPerLayer;
           n_pts = 100;
           graphzoom = 1;
           x = createMesh();
           h = self.computeHeights(x(:,1),x(:,2),n_pts,nF);
           figure(10)
           clf(10)     
           colorsc = ['r','g','b','c','m','y','k'];
           colorsc = fliplr(colorsc(1:self.data.nLabels));
           colorRGB = [1,0,0;0,1,0;0,0,1;0,1,1;1,0,1;1,1,0;0,0,0];
           colorRGB = flipud(colorRGB(1:self.data.nLabels,:));           
           switch type
               case 'contour'
                   for i = 1:nPL(end)
                       hold on
                       contour(x(:,1),x(:,2),h(:,:,i)',[0.5,0.5],'color',colorsc(i))
                   end
               case 'filled'
                   im = cell(size(self.data.Ytrain,2),1);
                   mymap = colormaps();
                   for i = 1:nPL(end)
                       hold on
                       ha(i) = axes;
                       im{i} = imagesc(x(:,1),x(:,2),h(:,:,end+1-i)');
                       im{i}.AlphaData = .3;                      
                       colormap(ha(i),mymap{i})
                       set(ha(i),'color','none','visible','off')
                       set(ha(i),'ydir','normal');
                   end
                   linkaxes(ha)
               case 'filledS'
                   [~,I] = max(h,[],3);
                   rgb = zeros(n_pts,n_pts,3);
                   for i = 1:n_pts
                       for j = 1:n_pts
                           for f = 1:nPL(end)                              
                               if I(i,j) == f
                                   rgb(j,i,:) = colorRGB(f,:);
                               end
                           end
                       end
                   end
                   RI = imref2d(size(I));
                   RI.XWorldLimits = [min(X(:,1)) max(X(:,1))];
                   RI.YWorldLimits = [min(X(:,2)) max(X(:,2))];
                   imshow(rgb,RI)
                   alpha(0.3);
           end  
           hold on
           title('Contour 0')
           self.data.plotdata(1,2);
           hold off

           function x = createMesh()
               gZ = graphzoom;
               extra_f1 = mean(X(:,1))*gZ;
               extra_f2 = mean(X(:,2))*gZ;
               x1 = linspace(min(X(:,1))-extra_f1,max(X(:,1))+extra_f1,n_pts)';
               x2 = linspace(min(X(:,2))-extra_f2,max(X(:,2))+extra_f2,n_pts)';
               x = [x1,x2];
           end

           function mymap = colormaps()
               colors = [1,0.3,0.3;   % r
                         0.3,1,0.3;   % g
                         0.3,0.3,1;   % b                         
                         0.3,1,1;     % c
                         1,0.3,1;     % m
                         1,1,0.3;     % y
                         0.3,0.3,0.3];% k
               mymap = cell(size(self.data.Ytrain,2),1);
               n = 100;
               w = ones(n/2,3);
               nLb = size(self.data.Ytrain,2);     
               for k = 1:nLb
                   mymap{k} = zeros(n,3);
                   r = linspace(1,colors(k,1),n/2)';
                   g = linspace(1,colors(k,2),n/2)';
                   bl = linspace(1,colors(k,3),n/2)';
                   mymap{k} = [w;r,g,bl];
               end
           end
        end

        function plotNetworkStatus(self)   
            layer = self.network.layer;
            nPL = self.neuronsPerLayer;
            nLy = length(nPL);
            neurons = cell(max(nPL),nLy);
            for i = 1:nLy-1
                if i == 1
                    maxTH = max(layer{i}.W(:));
                else
                    if maxTH < max(layer{i}.W(:))
                        maxTH = max(layer{i}.W(:));
                    end
                end
            end
            x_sep = 50;
            y_sep = 30;
            figure
            xlim([-20 nLy*x_sep-10])
            ylim([-20 max(nPL)*y_sep+10])
            set(gca,'XTick',[], 'YTick', [])
            box on
            hold on
            for i = 1:nLy
                for j = 1:nPL(i)
                    if i == 1
                        color = [1 0 0];
                    elseif i == nLy
                        color = [0 0.45 0.74];
                    else
                        color = [1 .67 .14];
                    end
                    prop.cx = (i-1)*x_sep;
                    prop.cy = (j-1)*y_sep;
                    prop.color = color;
                    prop.r = 10;
                    circlei = Circle(prop); %Circle composition
                    neurons{i,j} = circlei;
                    circlei.plot();
                end
            end
            r = [linspace(1,1,50)',linspace(0.75,0,50)',linspace(0.75,0,50)'];
            b = [linspace(0,0.75,50)',linspace(0,0.75,50)',linspace(1,1,50)'];
            rgb = [b;r];
            for i = 1:nLy-1
                for j = 1:nPL(i)
                    neuronb = neurons{i,j};
                    for k = 1:nPL(i+1)
                        neuronf = neurons{i+1,k};
                        wth = abs(layer{i}.W(j,k)/maxTH);
                        lw = 3*wth;
                        idx = round(wth*100);
                        if idx == 0
                            idx = 1;
                        end
                        linecolor = rgb(idx,:);
                        %linecolor = [0 , 0, 1];
                        line([neuronb.fx,neuronf.bx],[neuronb.fy,neuronf.by],'Color',linecolor,'LineWidth',lw)
                    end
                end
            end
            hold off
        end

        function drawConfusionMat(self)
            targets = self.data.Ytest;
            x = self.data.Xtest;
            outputs = self.network.getOutput(x);
            figure(1)
            plotconfusion(targets',outputs')
        end
    end

    methods (Access = private)
        function h_3D = computeHeights(self,x1,x2,n_pts,nF)
           nPL = self.neuronsPerLayer;
           X_test = zeros(n_pts,nF,n_pts);
           h = zeros(n_pts*nPL(end),n_pts);
           h_3D = zeros(n_pts,n_pts,nPL(end));
           for i = 1:n_pts
               x2_aux = ones(n_pts,1)*x2(i);
               xdata_test = [x1 , x2_aux];
               xful       = buildModel(xdata_test,self.data.polyGrade);
               X_test(:,:,i) = xful;
               h(:,i) = reshape(self.network.getOutput(X_test(:,:,i)),[n_pts*nPL(end),1]);
           end
           for j = 1:nPL(end)
               h_3D(:,:,j) = h((j-1)*n_pts+1:j*n_pts,:);
           end
       end 
    end
end