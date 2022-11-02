#{
----------------------- Cuadrados Minimos -------------------------
- Funcion de error: MSD_f(p) = 1/N.sum_i[0,n]((yi - f_p(xi))^2)
- Metodo derivando:
  1. Minimizamos por cada pj: -2/N.sum_i[0,n](yi - f_pj(xi)).(df/dpj)
  2. Armamos las matrices A y b derivando por cada parametro y resolvemos A.p = b
- Metodo automatico (algebraica):
  1. Queremos resolver A[#medidas,#parametros].p = b
  2. #medidas > #parametros => sobredeterminado => equivalente a minimizar abs(b - A.p)
  3. El problema es indeterminado, pero podemos usar el criterio de proyectar b en el espacio Im(A) = col(A)
  4. Para esto aprovechamos pinv(A) => p = pinv(A).b => calculamos pinv(A) mediante SVD
#}

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