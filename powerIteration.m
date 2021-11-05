1;
function [Xp, v] = powerIteration(A, Xi, error = 1e-5)
n = 1000;
Xi = Xi/norm(Xi);
vn = 1/eps;
  for i = 1:n
     Xp = A*Xi;
     v = Xp' * Xi;
     Xp = Xp/norm(A*Xi);
    if abs(v - vn) < error
      return;
    end
    vn = v;
    Xi = Xp;
  end
end

