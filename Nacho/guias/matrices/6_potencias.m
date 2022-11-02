#{
-------------------- Metodo de las Potencias -----------------------
- Busca el autovalor dominante (el de mayor modulo)
  1. v0 = (1, ..., 1)
  2. vi = A.versor(vi-1)
  3. li = vi.versor(vi-1)
- Criterios de convergencia: abs(A.versor(vi) - li.versor(vi)) < eps; abs(li - li-1) < eps;
- El autovalor mas chico, es el inverso del dominante de la inversa.
- Metodo inverso (sin calcular inversa):
  - El paso inductivo ahora es A.vi = versor(vi-1) y resolvemos el sistema
  - Para eso factoriamos A con PLU o QR, pues es solo una vez y amortiza
- Para sacar todos los autovalores podemos desfactorizarlos
- Desfactorizar autovalor i: B = A - li.versor(vi).versor(vi)'
#}

# -------------------- Potencias ----------------------- #
function [V, lambda, it, dif] = pow_met(A, v0, n=1000, tol=1e-10, tol_typ="lambda")

    V0 = v0;
    V0 = V0/norm(V0);

    V = [V0'];
    lambda = [];

    for it=1:n
        new_V = A*V(it,:)';
        new_lambda = dot(new_V, V(it,:));
        lambda = [lambda; new_lambda];
        new_V = new_V/norm(new_V);
        V = [V; new_V'];

        if it > 1
            if strcmp(tol_typ, "lambda") == 1
                dif = abs(new_lambda - lambda(it-1));
            elseif strcmp(tol_typ, "full") == 1
                dif = norm(A*new_V - lambda(it-1)*new_V);
            endif
            if dif < tol
                break
            endif
        endif
    endfor

endfunction

# -------------------- Potencias Inversas ----------------------- #
function [V, lambda, it, dif] = pow_met_inv(A, v0, n=1000, tol=1e-10, tol_typ="lambda")

    [Q, R] = qr(A);
    Q_inv = Q';

    V0 = v0;
    V0 = V0/norm(V0);

    V = [V0'];
    lambda = [];

    for it=1:n
        % R*new_V = Q_inv*V(it,:)

        new_V = upper_triang(R, Q_inv * V(it,:)');
        new_lambda = dot(new_V, V(it,:));
        lambda = [lambda; new_lambda];
        new_V = new_V/norm(new_V);
        V = [V; new_V'];

        if it > 1
            if strcmp(tol_typ, "lambda") == 1
                dif = abs(new_lambda - lambda(it-1));
            elseif strcmp(tol_typ, "full") == 1
                dif = norm(A*new_V - lambda(it-1)*new_V);
            endif
            if dif < tol
                break
            endif
        endif
    endfor

endfunction

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