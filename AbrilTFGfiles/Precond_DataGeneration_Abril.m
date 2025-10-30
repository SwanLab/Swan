
% DATA GENERATION

clc; clear; close all;

r = linspace(0,1,200); 
%r=0.4130;

K_all=zeros(8,8,length(r));
U_all=zeros(4901,19,length(r));

% Obtains the K coarse for each radius
guideElasticProblem_abril(0.1)

for j = 1:size(r,2)
    K = [];
    auxl = [];
    
    [~, u, l, mesh,Kcoarse] = LevelSetInclusionAuto_abril(r(j),1);

    K_all(:,:,j)=Kcoarse;

    % Reshapes U data and adds coordinates
    u_aux1=reshape(u.',8*2,[]).';
    u_aux2=[r(j)*ones(size(mesh.coord,1),1), mesh.coord, u_aux1];
    U_all(:,:,j)=u_aux2;

    %Designa un nom per cada linea corresponent a un radi
    string = strrep("UL_r"+num2str(r(j), '%.4f'), ".", "_")+"-20x20"+".mat"; 

    U         = u;
    L         = l;
    R         = r(j);
    K         = Kcoarse;

    % Guarda el workspace per cert radi
    %FileName=fullfile('AbrilTFGfiles','DataVariables',string)
    FileName=fullfile('AbrilTFGfiles','DataComparison',string);
    save(FileName, "U", "L", "K","mesh","R"); 
end



%% Reshapes the U data and saves it in a csv file

%reshaping coordinates since the Us are disposed as [ux1,uy1,ux2,uy2] for
%every coordinate
Udata=[];
for n=1:size(U_all,3)
    Udata=[Udata;U_all(:,:,n)];
end
T=array2table(Udata,"VariableNames",{'r','x','y','Tx1','Tx2','Tx3','Tx4','Tx5','Tx6','Tx7','Tx8' ...
    'Ty1','Ty2','Ty3','Ty4','Ty5','Ty6','Ty7','Ty8'});

uFileName = fullfile('AbrilTFGfiles', 'Udata.csv');
writematrix(Udata,uFileName);
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


