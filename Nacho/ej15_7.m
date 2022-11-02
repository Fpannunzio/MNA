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
theta=1;

%Dominio
xmin=0;
xmax=1;
Tmax=0.25;

%Resolucion
NIx=19;
Nt=1000;

%Condiciones de borde
bl=0;
br=0;

%Condicion iniciales
u0= @(x) 100*sin(pi*x);

%% alpha
a=1;

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
u=ones(Nt+1,NIx+2);
u(1,:) = u0(x);           %condicion inicial temperatura uniforme
u(:,1)=bl*ones(Nt+1,1);   %condicion de contorno
u(:,end)=br*ones(Nt+1,1); %condicion de contorno

%parametros del calculo

A = a*ht*theta/hx^2;
B = a*ht*(1-theta)/hx^2;

r1 = (1 + 2*A);
r2 = (1 -2*B);

%matriz n+1
d0np1 =r1*ones(NIx,1);
dsnp1 = -A*ones(NIx-1,1);
M=(diag(d0np1)+diag(dsnp1,-1)+diag(dsnp1,1));

%matriz n
d0n = r2*ones(NIx,1);
dsn = +B*ones(NIx-1,+1);
N=(diag(d0n)+diag(dsn,-1)+diag(dsn,1));


%   M*u^(n+1)=N*u^n +b; S= inv(M)*N
%     u^(n+1)=S*u^n + inv(M)*b
Minv= inv(M);
S   = Minv*N;
b=zeros(NIx,1);

%plot estado incial
plot(x,u(1,:))

%error en x=0.5
error=zeros(Nt, 1);

% sal == respuesta. sal(x, t)
sal=u(1,:)';
hold on
for k=2:Nt+1
    
    un=u(k-1,2:end-1)';
    
    b(1)    = A*u(k,1)      + B*u(k-1,1);
    b(NIx)  = A*u(k,NIx+2)  + B*u(k-1,NIx+2);
    
    %un+1
    u(k,2:end-1)=S*un+Minv*b;
    
    ux = u(k,floor((NIx+1)/2)+1);
    error(k-1) = abs(ux - 100*exp(-pi^2*k*ht)*sin(pi*0.5));
    
    if rem(k,np)==0
      plot(x,u(k,:))
      hold on
      sal=[sal u(k,:)'];
    endif
endfor

max(error)
min(error)
