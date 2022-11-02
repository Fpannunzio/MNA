%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Script para la resolución de la ecuacion de difusión 1D
                         u_tt=a^2 u_xx     
             con condiciones de borde períodicas 
                          u(xL,t)= gL(t) 
                          u(xR,t)= gR(t)
Se implementó para gR=br y gL=bl
             y condicion inicial
                          u(x,0)= cte3 
             solución analítica
                           u(x,t)=exp(-(x-ct)^2)            
      Se utiliza la discretización temporal 
           Euler progresivo  explicito theta=0
           Euler regresivo   implicito theta=1
           Crank Nicolson    semiimplicito  theta=1/2
           Rvieytes
           MNA 2C 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
clear all, clf

%Discrertizacion temporal
theta=0;

%Dominio
xmin=0;
xmax=10;
Tmax=80;

%Resolucion
NIx=4;
Nt=80;

%Condiciones de borde
bl=20;
br=100;

%Condicion iniciales
t0=30;

%% alpha cuadrado
a=2;

%% paso (se calcula solo)
hx=(xmax-xmin)/(NIx+1)
ht=Tmax/Nt

r=a*ht/hx^2

%cantidad de lineas ploteadas y soluciones exportadas
np=floor(Nt/5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Check estabilidad

if (theta==0.5) && r>1
    fprintf('Con los pasos temporal y espacial establecidos el esquema no converge\n')
    return
endif
if theta==0 && r>0.5
  fprintf('Euler progresivo no converge\n')
 % return
endif
%genera la grilla espacial y temporal
x=xmin:hx:xmax;
t=0:ht:Tmax;


%defino matriz de datos

%las temperaturas se almacenan en la matriz, cada tiempos es una fila.
%carga las condiones de contorno e iniciales
u=t0*ones(Nt+1,NIx+2);    %condicion inicial temperatura uniforme
u(:,1)=bl*ones(Nt+1,1);   %condicion de contorno izquierda
u(:,end)=br*ones(Nt+1,1); %condicion de contorno derecha

%parametros del calculo

% beta
A = r*theta;
% gamma
B = r*(1-theta);

r1 = (1 + 2*A);
r2 = (1 -2*B);

%matriz n+1
d0np1 = r1*ones(NIx,1);
dsnp1 = -A*ones(NIx-1,1);
M = diag(d0np1) + diag(dsnp1,-1) + diag(dsnp1,1);

%matriz n
d0n = r2*ones(NIx,1);
dsn = +B*ones(NIx-1,+1);
N=diag(d0n) + diag(dsn,-1) + diag(dsn,1);

% Teoria:
%   M*u^(n+1)=N*u^n +b; S= inv(M)*N
%     u^(n+1)=S*u^n + inv(M)*b
Minv= inv(M);
S   = Minv*N;
b=zeros(NIx,1);

%plot estado incial
plot(x,u(1,:))

% sal == respuesta. sal(x, t)
sal=u(1,:)';
hold on
for k=2:Nt+1
    
    %puntos interiores del estado anterior
    un=u(k-1,2:end-1)';
    
    b(1)=A*u(k,1)+B*u(k-1,1);
    b(NIx)=B*u(k-1,NIx+2)+A*u(k,NIx+2);
    
    %un+1
    u(k,2:end-1)=S*un+Minv*b;
    
    if rem(k,np)==0
      plot(x,u(k,:))
      hold on
      sal=[sal u(k,:)'];
    endif
endfor
