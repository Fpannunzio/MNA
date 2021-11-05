

w= 2*pi;
B= 0.5;
h= 0.01 %0.01 0.001

%Explicita
expl= h*[1/h 1; -w^2 (1/h)-B];

%Implicita
impl= (1/h)*inv([1/h -1;w^2 (1/h)+B]);

%Ratchet Y Clank
cn= inv([2/h -1;w^2 (2/h)+B]) * [2/h 1; -w^2 (2/h)-B];

A=cn;

y0= 1;
u0= 0;

tf= 2*pi;

v= [y0 u0]';

y= @(t) cos(w*t) * e^(t*B/2);

for t = 0:h:tf
  v=A*v;
end

df = v(1)
ex = y(tf)
