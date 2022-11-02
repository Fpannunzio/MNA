#{
----------------------------- Parabolica ---------------------------------
- Ecuacion del calor: Ut = a^2.Ux
- Drichlet: condiciones de contorno directas; Neumann: condiciones de contorno mediante derivadas
- r = a^2.ht/hx^2; b = o.r; l = (1-o).r
- Ecuacion general: -b.U_j-1_n+1 + (1+2b)U_j_n+1 - bU_j+1_n+1 = lU_j-1_n + (1-2l)U_j_n + lU_j+1_n  
- Metodos:
  1. FTCS: explicito - o = 1 - r < 1/2 (1/4 oscilaciones) - O(ht, hx^2)
  2. BTCS: implicito - o = 0 - incond. estable - O(ht, hx^2)
  3. CTCS: CN        - o = 1/2 - incond. estable - O(ht^2, hx^2)
- Condiciones de contorno Drichlet: U_j_0; U_0_n; U_N+1_n
  1. j = 1: (1+2b)U_1_n+1 - bU_2_n+1 = (1-2l)U_1_n + lU_2_n + (lU_0_n + bU_0_n+1)
  2. j = N: -b.U_N-1_n+1 + (1+2b)U_N_n+1 = lU_N-1_n + (1-2l)U_N_n + (lU_N+1_n + bU_N+1_n+1)
- Condiciones de contorno Neumann: U_j_0; Ux_0_n; Ux_N+1_n
  - Fantasmas
    1. U_-1_n = U_1_n - 2.hx.Ux_0_n
    2. U_N+2_n = U_N_n - 2.hx.Ux_N+1_n
  1. j = 0: (1+2b)U_0_n+1 - 2bU_1_n+1 = (1-2l)U_0_n + 2lU_1_n + 2hx(lUx_0_n + bUx_0_n+1)
  2. j = N+1: -2b.U_N_n+1 + (1+2b)U_N+1_n+1 = 2lU_N_n + (1-2l)U_N+1_n + 2hx(lUx_N+1_n + bUx_N+1_n+1)
- Solucion de la forma general M.U_n+1 = N.U_n + b
  - M asociado con b
  - N asociado con l
  - b termino independiente de los bordes
  - Sean A = inv(M).N y B = inv(M).b => U_n = A^n.U_j_0 + sum_i[0,n-1](A^i).B
#}

# ----------------------------- Parabolica del truchi ---------------------------------
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