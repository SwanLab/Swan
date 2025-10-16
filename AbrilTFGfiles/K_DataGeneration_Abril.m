
% DATA GENERATION

clc; clear; close all;

%r = linspace(0,1,200); 
r=0.4130;

K_all=zeros(8,8,length(r));

% Obtains the K coarse for each radius
guideElasticProblem_abril(0.1)

for j = 1:size(r,2)
    K = [];
    auxl = [];
    
    [~, u, l, mesh,Kcoarse] = LevelSetInclusionAuto_abril(r(j),1);
    K_all(:,:,j)=Kcoarse;

    %Designa un nom per cada linea corresponent a un radi
    string = strrep("UL_r"+num2str(r(j), '%.4f'), ".", "_")+"-20x20"+".mat"; 

    U         = u;
    L         = l;
    R         = r(j);
    K         = Kcoarse

    % Guarda el workspace per cert radi
    %FileName=fullfile('AbrilTFGfiles','DataVariables',string)
    FileName=fullfile('AbrilTFGfiles','DataComparison',string)
    save(FileName, "U", "L", "K","mesh","R"); 
end




%% Reshapes the data and saves it in a csv file

data=zeros(size(r,2),36);
for n=1:size(r,2)
    triangSup=triu(K_all(:,:,n));  %gets the triangular superior matrix
    clear row;
    row=[];
    for i=1:8
        for j=i:8
            row(end+1)=triangSup(i,j);
        end
    end
    data(n,:)=row;
end

data=[r.',data];

FileName = fullfile('AbrilTFGfiles', 'Kdata.csv');
writematrix(data,FileName);
