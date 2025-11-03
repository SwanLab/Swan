
% DATA GENERATION

clc; clear; close all;

r = 0:0.1:0.999; 
%r=0.4130;

K_all=zeros(8,8,length(r));
U_all1=zeros(761,19,length(r));
U_all2=zeros(761*8,5,length(r));

coarse_aux=repmat(eye(8),761*length(r),1);
%coarse_aux=join(compose('%d',coarse_aux),'');
%coarse=cell2mat(coarse_aux);
coarse=erase(string(num2str(coarse_aux)),' ');


% Obtains the K coarse for each radius
guideElasticProblem_abril(0.1)

for j = 1:size(r,2)
    
    [~, u, l, mesh,Kcoarse] = LevelSetInclusionAuto_abril(r(j),1);

    K_all(:,:,j)=Kcoarse;

    % Reshapes U data and adds coordinates
    u_aux1=reshape(u.',8*2,[]).';                                  %Joins the Tx and Ty coeff at the same line
    u_aux2=[r(j)*ones(size(mesh.coord,1),1), mesh.coord, u_aux1];  % Adds the radius and coordinates column
    U_all1(:,:,j)=u_aux2;                                          % Saves the result for each radius

    ux=reshape(u_aux1(:,1:8).',[],1);                       %Puts all the Tx coeff of each coord in a column concatenated
    uy=reshape(u_aux1(:,9:16).',[],1);                      %Idem with Ty
    coords=repelem(mesh.coord,8,1);                         %Coordinates column for each Tx and Ty 
    u_aux3=[r(j)*ones(size(coords,1),1),coords,ux,uy];    %concatenates the radius, coords and Tx Ty
    U_all2(:,:,j)=u_aux3;                                   %Saves the result
    
    
    %Designa un nom per cada linea corresponent a un radi
    string = strrep("UL_r"+num2str(r(j), '%.4f'), ".", "_")+"-20x20"+".mat"; 

    U         = u;
    L         = l;
    R         = r(j);
    K         = Kcoarse;

    % Guarda el workspace per cert radi
    FileName=fullfile('AbrilTFGfiles','DataVariables',string);
    %FileName=fullfile('AbrilTFGfiles','DataComparison',string);
    save(FileName, "U", "L", "K","mesh","R"); 
end



%% Reshapes the U data and saves it in a csv file

% Redimensioning the U_all1
Udata1=[];
for n=1:size(U_all1,3)
    Udata1=[Udata1;U_all1(:,:,n)];
end

% Redimensioning the U_all2
Udata2=[];
for n=1:size(U_all2,3)
    Udata2=[Udata2;(U_all2(:,:,n))];
end
Udata2=[coarse,Udata2];

T1=array2table(Udata1,"VariableNames",{'r','x','y','Tx1','Tx2','Tx3','Tx4','Tx5','Tx6','Tx7','Tx8' ...
    'Ty1','Ty2','Ty3','Ty4','Ty5','Ty6','Ty7','Ty8'});
uFileName = fullfile('AbrilTFGfiles', 'Udata1.csv');
writematrix(Udata1,uFileName);
%writetable(T1,uFileName);


T2=array2table(Udata2,"VariableNames",{'coarse','r','x','y','Tx','Ty'});
T2=T2(:,["r","x","y","coarse","Tx","Ty"]);

uFileName = fullfile('AbrilTFGfiles', 'Udata2.csv');
writematrix(Udata2,uFileName);
%writetable(T2,uFileName);


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


