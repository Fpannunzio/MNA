
conf = 2;

if conf == 1
  H  = [
    0   1/2 0     0;
    1/3 0   1/2   0;
    1/3 1/2 0     0;
    1/3 0   1/2   1;
  ];
  
elseif conf == 2
  H  = [
    0   1/2 0     0   1/5;
    1/3 0   1/2   0   1/5;
    1/3 1/2 0     0   1/5;
    1/3 0   1/2   0   1/5;
    0   0   0     1   1/5;
  ];
endif

n = 10;
b = 0.8;

N = size(H)(1);

% Hpr = b*H + (1-b)/N * ones(size(H));

Hpr = b*H + b/N * ones(N, 1) * [0 0 0 0 1] + (1-b)/N * ones(size(H));

pr = Hpr^n;

pr(:,1)
