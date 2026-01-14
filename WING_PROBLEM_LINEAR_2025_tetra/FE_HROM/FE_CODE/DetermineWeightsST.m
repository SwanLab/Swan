function [wST,wSTs] = DetermineWeightsST(detJeALL,weigREP,ngaus,g,nelem,wST,wSTs,nstrain,ROWSglo);
%%%%
% -------------------------------------------------------------------------
% COMMENTS (generated automatically by ChatGPT on 7-Nov-2025)
%
% PURPOSE:
%   Populate the vectors of quadrature weights scaled by Jacobian determinants.
%   Given detJ for all elements at Gauss point g and the Gauss weights,
%   this routine fills:
%     • wSTs : (nelem*ngaus x 1) with (detJ_e(g) * w_g) at the positions
%              corresponding to Gauss point g across all elements.
%     • wST  : (nelem*ngaus*nstrain x 1) with the same values as wSTs,
%              but each repeated nstrain times (block-wise), placed at the
%              indices specified by ROWSglo.
%
% INPUTS:
%   detJeALL : (nelem x 1) vector with determinant of J for Gauss point g
%              and all elements (e = 1..nelem).
%   weigREP  : (nelem x ngaus) matrix of Gauss weights replicated per element;
%              column j holds w_j for all elements.
%   ngaus    : number of Gauss points per element.
%   g        : current Gauss point index (1..ngaus).
%   nelem    : number of elements.
%   wST      : (nelem*ngaus*nstrain x 1) output vector to be filled (in/out).
%   wSTs     : (nelem*ngaus x 1) output vector to be filled (in/out).
%   nstrain  : number of strain components (e.g., 3 in 2D, 6 in 3D).
%   ROWSglo  : ((nelem*nstrain) x 1) linear indices selecting in wST the slots
%              that correspond to Gauss point g across all elements and all
%              nstrain repeats.
%
% OUTPUTS:
%   wSTs : updated with entries (detJ_e(g) * w_g) at indices g:ngaus:nelem*ngaus.
%   wST  : updated with nstrain-repeated copies of wSTs(g:ngaus:...) at ROWSglo.
%
% DETAILS / INDEXING:
%   - The pattern g:ngaus:nelem*ngaus selects, for the g-th Gauss point,
%     one position per element in the stacked (element-major) vector [e=1..nelem].
%     Concretely, positions are: g, g+ngaus, g+2*ngaus, ..., g+(nelem-1)*ngaus.
%   - wLOCa = detJeALL .* weigREP(:, g) computes the per-element product
%     detJ_e(g)*w_g. The code uses linear indexing weigREP(g:ngaus:nelem*ngaus),
%     which is equivalent to taking the g-th column across all rows.
%   - To fill wST (which stores nstrain copies for each element/Gauss pair),
%     the vector wLOCa (length nelem) is transposed, tiled nstrain times,
%     and reshaped column-wise into wLOCb of length (nelem*nstrain).
%     Those values are then written into wST at ROWSglo, which maps exactly
%     the slots of Gauss point g across all elements and their nstrain repeats.
%
% ASSUMPTIONS:
%   - detJeALL corresponds to the same Gauss point g provided as input.
%   - The calling code guarantees consistent stacking of elements/gauss points
%     in wSTs (element-major, gauss-minor) and a matching ROWSglo for wST.
%
% COMPLEXITY:
%   - Fully vectorized over elements; no per-element loops.
%
% SANITY:
%   - Negative or near-zero detJeALL values are not checked here; caller may
%     validate element orientations beforehand.
% -------------------------------------------------------------------------

wLOCa = detJeALL.*weigREP(g:ngaus:nelem*ngaus) ;
wSTs(g:ngaus:nelem*ngaus) =  wLOCa ;
%wLOCa = [5 6]' ;
wLOCb = repmat(wLOCa',nstrain,1) ;
wLOCb  = wLOCb(:) ;
wST(ROWSglo) = wLOCb(:);