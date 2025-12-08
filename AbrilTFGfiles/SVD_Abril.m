%SVD FOR THE NN OPTION 3

clc;
clear

r=0:0.02:0.999;

for i=1:size(r,2)
    string = strrep("UL_r"+num2str(r(i), '%.4f'), ".", "_")+"-50x50"+".mat"; 
    FileName=fullfile('AbrilTFGfiles','Data','50x50',string);
    
    load(FileName,"T","mesh");

    if i==1
        T_SVD=zeros(mesh.nnodes*mesh.ndim*8,size(r,2));
    end

    T_SVD(:,i)=T(:);

end

[U,S,V]=svd(T_SVD,'econ');

%% Graphics


% PLot the V columns grouped in 10
step=10;
Nwindow=ceil(size(V,2)/step);
idx=1;

for j=1:Nwindow
    figure('Position',[75 100 1400 600]);
    tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
    for i=1:step
        ax=nexttile;
        plot(r,V(:,idx), 'LineWidth', 1.5);
        xlabel('r');
        ylabel("V(:,"+idx);
        title("V-"+ idx+" ; r="+r(1,idx));
        grid on
        idx=idx+1;
    end
end

tiledlayout(10,5,'TileSpacing','compact','Padding','compact');
for i=1:50
    ax=nexttile;
    plot(r,V(:,i), 'LineWidth', 1.5);
    xlabel('r');
    ylabel("V"+i);
    title("V"+ i+" ; r="+r(1,i));
    grid on
    i=i+1;
end

% Plot the S
figure
plot(log(diag(S)),'LineWidth',1.5);
title("S singular values");
ylabel('value');
xlabel('column');

% Plot U
%s.mesh=mesh;
%s.order='P1';
%n=10;
%
%Ufun=cell(1,n);
%for i=1:n
%    s.fValues=U(:,i);
%    Ufun{1,i}=LagrangianFunction(s);
%    Ufun{1,i}.plot();
%end

figure('Position',[75 100 1400 600]);
tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
for i=1:10
    ax=nexttile;
    plot(U(:,i));
    ylabel("U values");
    title("U column "+ i);
end


    

