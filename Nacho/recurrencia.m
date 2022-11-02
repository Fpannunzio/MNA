
conf = 5;

% Ej 8,1
if conf == 1
  M = [
    0.8 0.3;
    0.2 0.7
  ];

  x0 = [0 5]';

% Ej 8,2  
elseif conf == 2
  M = [
    0.5     2/5;
    -13/125 11/10
  ];

  x0 = [0 5]';
  
% Ej 8,3
elseif conf == 3
  M = [
    0.4  0.5;
    0.2  0.1
  ];

  x0 = [0 5]';
  
  % El crecimiento debe ser que el autovalor mas grande da 1.02
  % El ratio sale de hacer M(1,:) / M(2,:)
  
% Ej 8,4
elseif conf == 4
  M = [
    0.5  2/3;
    0.5  1/3
  ];

  x0 = [1 0]';
  
% Ej 8,5
elseif conf == 5
  M = [
    0   0   0   0   1/2   0;
    1/4 0   0   0   0     0;
    0   1/2 0   0   0     0;
    1/4 1/2 0   0   1/2   0;
    1/4 0   1   1   0     1;
    1/4 0   0   0   0     0;
  ];

  x0 = [1 0]';
  
endif
n = 100;

[P, D] = eig(M);

Mn = P * D.^n * inv(P);

real(Mn)

%x = Mn * x0