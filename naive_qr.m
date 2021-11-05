1;
function [Q,R] = naive_qr(A)

    R = A;
    Q = 0;

    for i=1:columns(A)

        x = R(:,1);
        s = sign(x(1))*norm(x);

        % aux = zeros(ndims(x),1);
        % aux(1) = 1;

        % u = x + s*aux

        u = x;
        u(1) += s;

        H = eye(rows(u), rows(u)) - (2/norm(u)**2) * u*transpose(u);

        if i == 1
            Q = H;
        else
            BASE = eye(columns(Q), columns(Q));

            for r=i:rows(Q)
                for l=i:columns(Q)
                    BASE(r,l) = H(r-i+1,l-i+1);
                endfor
            endfor

            Q = Q*BASE;
        endif

        if i < columns(A)
            HA = H*R;
            R = HA(2:end,2:end);
        endif

    endfor

    Q;
    R = transpose(Q)*A;

endfunction