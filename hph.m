% C = hph(H,P)
% 
% Gives the first column (rest redundant) of H*P*H'.
% Testing for N = 2000 and M = 200, This method is ~12 times faster than
% multiplying the matrices (doing H*P*H') and extracting the first column.
% 
% IN:   H = NxM model ("regressor") matrix
% IN:   P = MxM prior, must be diagonal
% 
% OUT:  C = Nx1, the first column of H*P*H'
% 
function C = hph(H,P)
HH = H.*H(1,:);
C = HH*diag(P);
end
