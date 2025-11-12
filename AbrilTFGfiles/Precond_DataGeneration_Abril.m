
% DATA GENERATION

clc; clear; close all;

r = 0:0.1:0.999; 
%r=0.5;


T_all=zeros(761,19,length(r));

doplot=false();

for j = 1:size(r,2)
    
    [~, u, l, mesh,Kcoarse] = LevelSetInclusionAuto_abril(r(j),1,doplot);

    % Initialization for K_all and T_all
    if j==1
        K_all=zeros(8,8,length(r));
        T_all=zeros(mesh.nnodes,19,length(r));
    end

    K_all(:,:,j)=Kcoarse;
    
    % Reshapes U data and adds coordinates
    t1=reshape(u(:,1).',2,[]).';     % Joins the Tx and Ty coeff at the same line
    t2=reshape(u(:,2).',2,[]).';                                  
    t3=reshape(u(:,3).',2,[]).';                                  
    t4=reshape(u(:,4).',2,[]).';                                  
    t5=reshape(u(:,5).',2,[]).';                                  
    t6=reshape(u(:,6).',2,[]).';   
    t7=reshape(u(:,7).',2,[]).';
    t8=reshape(u(:,8).',2,[]).';

    t_aux=[r(j)*ones(size(mesh.coord,1),1), mesh.coord, t1,t2,t3,t4,t5,t6,t7,t8];  % Adds the radius and coordinates column

    T_all(:,:,j)=t_aux;   % Saves the result for each radius


    %Designa un nom per cada linea corresponent a un radi
    string = strrep("UL_r"+num2str(r(j), '%.4f'), ".", "_")+"-20x20"+".mat"; 

    T         = u;
    L         = l;
    R         = r(j);
    K         = Kcoarse;

    % Guarda el workspace per cert radi
    FileName=fullfile('AbrilTFGfiles','DataVariables',string);
    %FileName=fullfile('AbrilTFGfiles','DataComparison',string);
    save(FileName, "T", "L", "K","mesh","R"); 
end



%% Reshapes the U data and saves it in a csv file

% Redimensioning the U_all1
TData=[];
for n=1:size(T_all,3)
    TData=[TData;T_all(:,:,n)];
end

T=array2table(TData,"VariableNames",{'r','x','y','Tx1','Ty1','Tx2','Ty2','Tx3','Ty3','Tx4','Ty4' ...
    'Tx5','Ty5','Tx6','Ty6','Tx7','Ty7','Tx8','Ty8'});

uFileName = fullfile('AbrilTFGfiles', 'Tdata.csv');
writematrix(TData,uFileName);
%writetable(T,uFileName);



%% Reshapes the K data and saves it in a csv file

kdata=zeros(size(r,2),36);
for n=1:size(r,2)
    triangSup=triu(K_all(:,:,n));  %gets the triangular superior matrix
    clear row;
    row=[];
    for i=1:8
        for j=i:8
            row(end+1)=triangSup(i,j);
        end
    end
    kdata(n,:)=row;
end

kdata=[r.',kdata];
kFileName = fullfile('AbrilTFGfiles', 'Kdata.csv');
writematrix(kdata,kFileName);


