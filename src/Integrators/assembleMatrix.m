function A = assembleMatrix(Aelem,f1,f2)
dofsF1 = f1.getDofConnec();
if isequal(f1, f2)
    dofsF2 = dofsF1;
else
    dofsF2 = f2.getDofConnec();
end
nDofs1     = numel(f1.fValues);
nDofs2     = numel(f2.fValues);
ndofsElem1 = size(Aelem, 1);
ndofsElem2 = size(Aelem, 2);

[iElem, jElem] = meshgrid(1:ndofsElem1, 1:ndofsElem2);
iElem = iElem(:);
jElem = jElem(:);

dofsI = dofsF1(:, iElem);
dofsJ = dofsF2(:, jElem);

rowIdx = dofsI(:);
colIdx = dofsJ(:);
Aval   = permute(Aelem,[3 2 1]);
values = Aval(:);
A = sparse(rowIdx, colIdx, values, nDofs1, nDofs2);
end


