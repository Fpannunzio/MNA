1;
function [Xp, v] = powerIterationTryingValues(A, error=1e-5)
n = 50;
Xp(:, 1) = A(:,1);
v(1) = 0;
for i = 1:n
  Xi = (2^(-i))*[1 -3 1]' + [1 0 1]';
  Xi = Xi/norm(Xi);
  vn = 1/eps;
    for k = 1:n*20
       Xp(:,i) = A*Xi;
       v(i) = Xp(:,i)' * Xi;
       Xp(:,i) = Xp(:,i)/norm(A*Xi);
      if abs(v(i) - vn) < error
        break;
      end
      vn = v(i);
      Xi = Xp(:,i);
    end
  end
end