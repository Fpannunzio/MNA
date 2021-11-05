1;
function x = lower_triang(L, b)

    x = zeros(columns(L), 1);
    n = rows(L);

    x(1) = b(1)/L(1, 1);

    for idx=2:n
        s = 0.0;
        for jota=1:idx-1
            s = s + L(idx, jota)*x(jota);
        endfor
        x(idx) = (b(idx) - s)/L(idx, idx);
    endfor

    % return x

endfunction