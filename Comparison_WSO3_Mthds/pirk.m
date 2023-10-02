function [A, b, c] = pirk(cTilde)
    p = length(cTilde) - 1;
    V = fliplr(vander(cTilde));
    S = diag(1 ./ (1:p), -1);
    
    svInv = S / V;
    ATilde = V * svInv;
    %bTilde = svInv(end, :).';
    bTilde = sum(svInv)';
    A = kron(diag(ones(p - 1, 1), -1), ATilde);
    b = [zeros(p^2 - 1, 1); bTilde];
    b = b';
    c = sum(A, 2);
end