#{
----------------------------- Ejemplo IVP --------------------------------
- f' = -a.y; y(0) = y0
  1. Elegimos Aprox y recursion: ej. forward/explicito: (yn+1 - yn)/h = -a.yn
  2. yn+1 = y0.(1 - a.h)^n
  3. Confimrar estabilidad: abs(yn+1/yn) <= 1

----------------------------- Ejemplo BVP --------------------------------
- y' + 2.x.y = 0
- Centrado-explicito: (yn+1 - yn-1)/(2.h) + 2.x.yn = 0
  - -yn-1 + 4.h.x.yn + yn+1 = 0
  - Armas la matriz evaluando con todos los n
#}

# ----------------------------- IVP -------------------------------- #
{
  # ---------------------- PARAMETROS ------------------------------
  sigma   = 1/2;                    # 1/2 = centrado, 1 = implicita, 0 = explicita
  h       = 0.1;                    # El paso
  t_init  = 0;                      # cota izquierda del intervalo
  t_end   = 0.4;                    # cota derecha del intervalo
  N       = (t_end - t_init) / h;   # cantidad de puntos
  y0      = 1;                      # condicion inicial, y(t_init) = y0

  next_y = @(y, t) (h*t**2 + (1 - h*(1 - sigma)) * y) / (1 + h*sigma); # Usando forward

  # -------------------------- Iteramos -----------------------------------

  y = [y0];
  for i = 2:N + 1
    y(i) = next_y(y(i - 1), t_init + i*h);
  endfor

  y
}

# ----------------------------- BVP -------------------------------- #
{
  alpha = 2;
  h = 0.05;

  a = -(2 + alpha**2*h**2);

  M = diag([1, 1, 1, 2], -1) + diag(a*ones(5, 1)) + diag(ones(4, 1), 1)
  b = [-1, 0, 0, 0, -2*h];

  b*inv(M)
  cosh(alpha * (1 - (h*[1:5]))) / cosh(alpha)
}

# ----------------------------- Padding ---------------------------------
{
  %{
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              Script para la resolución de la ecuacion de difusión 1D
                           u_tt=a^2 u_xx     
               con condiciones de borde períodicas 
                            u'(xL,t)= gL(t) 
                            u'(xR,t)= gR(t)
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
  theta=1/2;

  %Dominio
  xmin=0;
  xmax=10;
  Tmax=100;

  %Resolucion
  NIx=4;
  Nt=100;

  %Condiciones de borde
  ul=20
  ur=100
  m=80/10
  bl=m;
  br=m;
  gl= @(t) bl
  gr= @(t) br

  %Condicion iniciales
  t0=30;

  %% alpha
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

  NIx = NIx + 2; %ajustado para Newman

  %defino matriz de datos

  %las temperaturas se almacenan en la matriz, cada tiempos es una fila.
  %carga las condiones de contorno e iniciales

  u=t0*ones(Nt+1,NIx);    %condicion inicial temperatura uniforme
                          % con Newman solo hay puntos interiores
  %Solo en t0 EN NEWMAN
  u(1,1)    =ul;   %condicion de contorno
  u(1,end)  =ur;   %condicion de contorno

  %parametros del calculo

  A = a*ht*theta/hx^2;
  B = a*ht*(1-theta)/hx^2;

  r1 = (1 + 2*A);
  r2 = (1 -2*B);

  %matriz n+1
  d0np1 =r1*ones(NIx,1);
  dsnp1 = -A*ones(NIx-1,1);
  M=(diag(d0np1)+diag(dsnp1,-1)+diag(dsnp1,1));

  %ajuste para Newman
  M(1,2)        = -2*A;
  M(NIx,NIx-1)  = -2*A;

  %matriz n
  d0n = r2*ones(NIx,1);
  dsn = +B*ones(NIx-1,+1);
  N=(diag(d0n)+diag(dsn,-1)+diag(dsn,1));

  %ajuste para Newman
  N(NIx,NIx-1)  = 2*B;
  N(1,2)        = 2*B;


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
      
      un=u(k-1,:)';  
      
      b(1)    = -2*hx*(A*gl(k)  + B*gl(k-1));
      b(NIx)  =  2*hx*(A*gr(k)  + B*gr(k-1));
      
      %un+1
      u(k,:)=S*un+Minv*b;
      
      if rem(k,np)==0
        plot(x,u(k,:))
        hold on
        sal=[sal u(k,:)'];
      endif
  endfor
}