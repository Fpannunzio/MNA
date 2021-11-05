function x = upper_triang(U, b)

    x = zeros(columns(U), 1);
    n = rows(U);

    x(n) = b(n)/U(n, n);

    for idx=n-1:-1:1
        s = 0.0;
        for jota=idx+1:n
            s = s + U(idx, jota)*x(jota);
        endfor
        x(idx) = (b(idx) - s)/U(idx, idx);
    endfor

    % return x

endfunction