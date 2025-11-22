%SVD FOR THE NN OPTION 3

clc;
clear

r=0:0.05:0.999;

for i=1:size(r,2)
    string = strrep("UL_r"+num2str(r(i), '%.4f'), ".", "_")+"-20x20"+".mat"; 
    FileName=fullfile('AbrilTFGfiles','DataVariables','20x20',string);
    
    load(FileName,"T","mesh");

    if i==1
        T_SVD=zeros(mesh.nnodes*mesh.ndim*8,size(r,2));
    end

    T_SVD(:,i)=T(:);

end

[U,S,V]=svd(T_SVD);

