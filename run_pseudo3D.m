% along-strike pressure diffusion with "local" elasticity
% and rate-and-state friction
%
% Eric M. Dunham

% load problem set up parameters (uncomment)

%params0 % nominal parameters
%params1 % higher perm
%params2 % lower perm
params3 % nominal parameters, stop injecting after 3 yr but simulate 5 yr

% discretization
Dxx = secondDerivative(M.nx,M.h);
delta = zeros(M.nx-1,1); delta(M.i-1) = 1/M.h;

% ta = time vector (a=adaptive)
% Da,Psia,Va,taua,pa = solution (a=adaptive)

t=0;

% store solution
ta=t; Da=D0; Psia=Psi0; pa=p0;
[V1,G1,tau1, Neff1] = sliderODEvector(Da(:,end),Psia(:,end),pa(:,end),t,M);
Va=V1; taua=tau1; Neffa=Neff1; % stage 1 values are stored

err=0; dta=dt;

while t<tmax
    
    % adjust dt to stop at tmax
    if t+dt>tmax, dt=tmax-t; end
    
    % three-stage method with embedded error estimate
    
    [V2,G2] = sliderODEvector(Da(:,end)+0.5*dt*V1,Psia(:,end)+0.5*dt*G1,pa(:,end),t+0.5*dt,M);
    [V3,G3] = sliderODEvector(Da(:,end)+dt*(-V1+2*V2),Psia(:,end)+dt*(-G1+2*G2),pa(:,end),t+dt,M);
    
    % second order update
    D2   = Da(:,end)  +dt/2*(V1+V3);
    Psi2 = Psia(:,end)+dt/2*(G1+G3);
    
    % third order update
    D3   = Da(:,end)  +dt/6*(V1+4*V2+V3);
    Psi3 = Psia(:,end)+dt/6*(G1+4*G2+G3);
    
    q = 2; % order of accuracy of lower order update
    
    % local error estimate
    er = norm([D2-D3; Psi2-Psi3]);

    if er<tol
        % update solution
        t = t+dt; ta=[ta; t];
        Da = [Da D3]; Psia = [Psia Psi3]; % use third-order update
        
        % store error and time step
        err=[err; er]; dta=[dta; dt];

        % update pressure (diffusion solve) and store
        A = speye(M.nx-1,M.nx-1)-dt*M.c*Dxx;
        s = dt/(M.beta*M.phi)*delta;
        s(M.i-1) = s(M.i-1).*M.q0'.*heaviside(t-M.ti)'.*heaviside(M.tstop-t)';
        p = [0; A\(pa(2:M.nx,end)+s); 0];
        pa = [pa p];
        
        % evaluate stage 1 values for next time step
        [V1,G1,tau1,Neff1] = sliderODEvector(Da(:,end),Psia(:,end),pa(:,end),t,M);
        Va=[Va V1]; taua=[taua tau1]; Neffa=[Neffa Neff1]; % stage 1 values are stored
        
        if dispAcceptReject, disp(['accept ' num2str(t) ' ' num2str(er) ' ' num2str(dt)]), end
        if plotSolDuringSim
            figure(1),clf
            subplot(3,1,1)
            semilogy(M.x*1e-3,Va(:,end))
            subplot(3,1,2)
            plot(M.x*1e-3,pa(:,end)*1e-6)
            subplot(3,1,3)
            plot(M.x*1e-3,Da(:,end))
            drawnow
        end
        
    else
        if dispAcceptReject, disp(['reject ' num2str(t) ' ' num2str(er) ' ' num2str(dt)]), end
    end
    
    % adjust time step
    dt = safety*dt*(tol/er)^(1/(q+1));
    dt = min(dt,dtmax);

end

% plot
figure(1),clf
subplot(3,1,1)
semilogy(ta/oneyear,Va(M.i,:))
xlabel('time (yr)')
ylabel('V (m/s)')

subplot(3,1,2)
plot(ta/oneyear,taua(M.i,:),ta/oneyear,Neffa(M.i,:))
xlabel('time (yr)')
ylabel('MPa')
legend('\tau','\sigma-p')

subplot(3,1,3)
plot(ta/oneyear,Da(M.i,:))
xlabel('time (yr)')
ylabel('slip (m)')

figure(2)
semilogx([1e-14 1e-4],M.f0+(M.a-M.b)*log([1e-14 1e-4]/M.V0))
hold on
for i=1:length(M.i)
    semilogx(Va(M.i(i),:),taua(M.i(i),:)./Neffa(M.i(i),:))
end
hold off
xlabel('V (m/s)')
legend('f_{ss}','f','location','southeast')

figure(3)
contourInterval = 0.5*oneyear; tInterval = [0:contourInterval:tmax];
nInterval = length(tInterval); Dplot = nan(M.nx+1,nInterval); pplot = nan(M.nx+1,nInterval);
for i=1:M.nx+1
    Dplot(i,:) = interp1(ta',Da(i,:),tInterval);
    pplot(i,:) = interp1(ta',pa(i,:),tInterval);
end
subplot(2,1,1)
plot(M.x*1e-3,Dplot,'b')
xlabel('x (km)')
ylabel('slip (m)')
subplot(2,1,2)
plot(M.x*1e-3,pplot*1e-6,'b')
xlabel('x (km)')
ylabel('pressure change (MPa)')

% functions below here

function [V,G,tau,Neff] = sliderODEvector(D,Psi,p,t,M)

    V = nan(M.nx+1,1);
    G = nan(M.nx+1,1);
    tau = nan(M.nx+1,1);
    Neff = nan(M.nx+1,1);
    
    for i=1:M.nx+1
        [V(i),G(i),tau(i),Neff(i)] = sliderODE(D(i),Psi(i),p(i),t,M.x(i),M);
    end
    
end

function [V,G,tau,Neff] = sliderODE(D,Psi,p,t,x,M)
    
    % evaluate stress when V=0

    tauLock = M.tau0+M.dtaudt*t-M.K*D;

    % set bounds on V for root-finding

    if tauLock>0
        Vmin = 0; Vmax = tauLock/M.eta;
    else
        Vmin = tauLock/M.eta; Vmax = 0;
    end

    % evaluate effective normal stress

    Neff = M.N-M.dpdt*t-p*1e-6;
    if any(Neff<=0), disp('tensile effective normal stress'), end
    
    % solve stress=strength for V
    
    atol = 1e-14; rtol = 1e-6;
    
    V = hybrid(@(V) solveV(V,tauLock,Neff,Psi,M) ,Vmin,Vmax,atol,rtol);
    
    % then evaluate tau

    tau = tauLock-M.eta*V;

    % and state evolution, G = dPsi/dt

    if V==0
        G = 0; % special case to avoid log(0)
    else
        f = tau/Neff; fss = M.f0+(M.a-M.b)*log(V/M.V0);
        G = -V/M.dc*(f-fss);
    end
        
end


function residual = solveV(V,tauLock,Neff,Psi,M)

    stress = tauLock-M.eta*V;
    f = M.a*asinh(V/(2*M.V0)*exp(Psi/M.a));
    strength = f*Neff;
    residual = stress-strength;
    
end


function [x,err]=hybrid(func,a,b,atol,rtol)

  % hybrid method solves func(x)=0 for some root x within (a,b)
  % returns x, estimate of root with absolute error less than atol
  % or relative error less than rtol
  
  % function values at endpoints
  fa = func(a);
  fb = func(b);

  % make sure root is bracketed; otherwise return
  if sign(fa)==sign(fb) | isnan(fa) | isnan(fb)
    disp('error: root not bracketed or function is NaN at endpoint')
    x = NaN; err = NaN;
    return
  end

  % set up secant method, storing old values as xold and fold, new ones as x and f
  % use bisection brackets to start secant (this is somewhat arbitrary)
  xold = a;
  fold = fa;
  x = b;
  f = fb;
  
  % begin iterations,
  % keeping track of error at each iteration in vector err
  n = 0; err = [];
  update = 'input'; % character string stating type of update used in previous interation
  while b-a>atol+rtol*abs(x) % safe to have infinite loop since bisection guaranteed to converge
      
      err = [err b-a]; % add to end of vector the current error (interval width)

      % formatted printing so you can watch method converge
      %fprintf('%6i %20.10f %20.10f %20.10f %s\n',n,a,x,b,update)

      n = n+1; % iteration number

      % first calculate (tenative) secant update
      dfdx = (f-fold)/(x-xold); % approximation to df/dx
      dx = -f/dfdx; % update to x
      xs = x+dx; % secant update
      
      % determine if secant method will be used
      if (xs<a) | (xs>b)
          use_secant = false;  % not if update outside (a,b)
      else
          fs = func(xs); % function value at secant update
          % calculate interval reduction factor = (old interval width)/(new interval width)
          if sign(fs)==sign(fa)
              IRF = (b-a)/(b-xs); % would update a=xs
          else
              IRF = (b-a)/(xs-a); % would update b=xs
          end
          if IRF<2
              use_secant = false;
          else
              use_secant = true;
          end
      end

      xold = x; fold = f; % store these values for next iteration

      % now update
      if use_secant
          update = 'secant';
          x = xs;
          f = fs;
      else
          update = 'bisection';
          x = (a+b)/2; % midpoint
          f = func(x); % function value at midpoint
      end
      
      % update one endpoint based on sign of function value at updated x
      if sign(f)==sign(fa)
          a = x;
      else
          b = x;
      end
      
  end
  
end

function D2 = secondDerivative(N,h)

    % second derivative matrix:

    % first create array holding nonzero diagonals
        
    % (N-1) by 3 array of coefficients
    e = ones(N-1,1);
    a = [e -2*e e]/h^2;
    
    % two entries are not used (because super- and 
    % sub-diagonals only have N-2 entries, not N-1)
    a(end,1) = NaN; % not used, can be set to any value
    a(1  ,3) = NaN; % not used, can be set to any value
    
    D2 = spdiags(a,[-1 0 1],N-1,N-1); % sparse (N-1) by (N-1) array

end

