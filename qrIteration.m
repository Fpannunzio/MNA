1;
function [Xp, V] = qrIteration(A, n=100)
  Xp = eye(size(A)(1));
  V = A;
  for i = 1:n
    [Q R] = qr(V)
    Xp = Xp*Q;
    V = R*Q
  endfor
end