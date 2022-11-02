#{
------------- Reflector Housholder ---------------
- H_l(v) = v - 2.<l.v>.l con l versor (sino dividido <l.l>)
- H_l = I - 2.l.l' con l versor
- H_l es simetrica y ortonormal
- Reflejar x solo a la componente i:
  1. s = sg(xi)*norm(xi)
  2. u = x + s.ei (le sumamos s en la componente i)
  3. H_u(x) = a.ei (vector solo con componente i)

-------------------- Factoriacion QR -----------------------
- A = QR (A[m,n])
- Q ortonormal por columnas
- R triangular superior

- Version Completa: Q[m,m], R[m,n]
  1. Ortonormalizamos las columnas de A con Grand Schmidt para obtener Q
  2. Como inv(Q) = Q' => R = Q'A

- Version Reducida (housholder): Q[m,n], R[n,n]
  1. R1 = A; H1 = H_u de reflexion de la primera componente de R1
  2. Ri = (Hi-1 * Ri-1)[2:end, 2:end] (sacando la primer col/row); Hi = H_u de reflexion de la primera componente de Ri
  3. Repetir hasta que Ri tenga dimension 1
  4. Q = sum(Hi) rellenandolos con 0 a izquierda
  5. Como inv(Q) = Q' => R = Q'A
#}

# ------------- QR Reducido (Housholder) --------------- #
function [Q,R] = naive_qr(A)

    R = A;
    Q = 0;

    for i=1:columns(A)

        x = R(:,1);
        s = sign(x(1))*norm(x);

        % aux = zeros(ndims(x),1);
        % aux(1) = 1;

        % u = x + s*aux

        u = x;
        u(1) += s;

        H = eye(rows(u), rows(u)) - (2/norm(u)**2) * u*transpose(u);

        if i == 1
            Q = H;
        else
            BASE = eye(columns(Q), columns(Q));

            for r=i:rows(Q)
                for l=i:columns(Q)
                    BASE(r,l) = H(r-i+1,l-i+1);
                endfor
            endfor

            Q = Q*BASE;
        endif

        if i < columns(A)
            HA = H*R;
            R = HA(2:end,2:end);
        endif

    endfor

    Q;
    R = transpose(Q)*A;

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