1;
function [v, eig] = powerIteration(A, v, error = 1e-5)
  n = 1000;
  prev_v = v/norm(v);
  prev_eig = 1/eps;
  for i = 1:n
     v = A*prev_v;
     eig = v' * prev_v;
     v = v/norm(v);
    if abs(eig - prev_eig) < error
      return;
    end
    prev_eig = eig;
    prev_v = v;
  end
end

