1;
function [Xp, v] = exercise9_6(error = 1e-6)
n = 15
v(1) = 0
for i = 1:n
  t=8
  A = [0 -10*i; 1 (10 + i)];
  Xi = [1 1]'
  Xi = Xi/norm(Xi)
  vn = 1/eps
    for k = 1:n*200
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