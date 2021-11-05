1;
function [V, lambda, it, dif] = pow_met(A, v0, n=1000, tol=1e-10, tol_typ="lambda")

    V0 = v0;
    V0 = V0/norm(V0);

    V = [V0'];
    lambda = [];

    for it=1:n
        new_V = A*V(it,:)';
        new_lambda = dot(new_V, V(it,:));
        lambda = [lambda; new_lambda];
        new_V = new_V/norm(new_V);
        V = [V; new_V'];

        if it > 1
            if strcmp(tol_typ, "lambda") == 1
                dif = abs(new_lambda - lambda(it-1));
            elseif strcmp(tol_typ, "full") == 1
                dif = norm(A*new_V - lambda(it-1)*new_V);
            endif
            if dif < tol
                break
            endif
        endif
    endfor

endfunction