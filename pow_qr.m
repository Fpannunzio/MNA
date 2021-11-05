1;
function [L, k, dif] = pow_qr(A, numiter=1000, tol=1e-10)

    M = A;
    L = [diag(A)'];
    for k=1:numiter
        [Q R] = qr(M);
        M = R*Q;
        D = diag(M)';
        L = [L; D];
        dif = sum(abs(L(k,:) - D));
        if dif < tol
            break
        endif
    endfor

endfunction