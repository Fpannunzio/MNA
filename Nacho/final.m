%Dominio
xmin=0;
xmax=1;

y0 = 0
y1=0

%Resolucion
NIx=20

D=0.1;

h=(xmax-xmin)/(NIx+1);


d0np1 = 4*D*ones(NIx,1);
dsnp1 = (h-2*D)*ones(NIx-1,1);
dsnp2 = (-h-2*D)*ones(NIx-1,1);
M = diag(d0np1) + diag(dsnp2,-1) + diag(dsnp1,1);

b=2*h*h*ones(NIx,1);
b(1)=2*h*h + (h+2*D)*y0;
b(end)=2*h*h - (h-2*D)*y1;

[Q, R] = qr(M);

b=Q'*b;

x = zeros(columns(R), 1);
n = rows(R);

x(n) = b(n)/R(n, n);

for idx=n-1:-1:1
    s = 0.0;
    for jota=idx+1:n
        s = s + R(idx, jota)*x(jota);
    endfor
    x(idx) = (b(idx) - s)/R(idx, idx);
endfor




f=@(t) t - (exp(-(1-t)/D)-exp(-1/D))/(1-exp(-1/D));

sum((f(0:h:1) - [0 x' 0]).^2)/(NIx+2)

