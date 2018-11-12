function openmatlab(nLAB)

nLabImp_TRAINING = nLAB;

nLab = matlabpool('size');
if nLab>0
    disp(['Se utiliza la configuracion activa: ',get(findResource(),'configuration'),...
        ' de ',num2str(nLab),' laboratorios.']);
else
    %Recordar que para pocos elementos no conviene paralelizar.
    if nLabImp_TRAINING==0
        matlabpool open
    elseif nLabImp_TRAINING>0
        matlabpool('open',nLabImp_TRAINING)
    end
end


end
