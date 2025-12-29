function T_trained=computeT_NN(mesh,R,T_NN,pol_deg)
    T_trained=[];
        for i=1:size(mesh.coord,1)  % Evaluates all the coordenates
            dataInput=[R,mesh.coord(i,:)];
            dataFull=Data.buildModel(dataInput,pol_deg);
            Taux1=T_NN.computeOutputValues(dataFull).';
            Taux2=reshape(Taux1,2,[]);
            T_trained=[T_trained;Taux2];
        end
end