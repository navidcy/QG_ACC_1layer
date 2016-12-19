clearvars ?global;
clc;

global eta etahat detadx detady del2 ksq KX KY Nx Ny F beta mu nu4 nu2 FILTER;

% grid points
Nx=256;
Ny=Nx;

% domain size
Lx=2*pi;Ly=2*pi;

% linear drag coefficient
mu  = 1e-2;

% diffusion
nu2  = 0;

% hyper-diffusion
nu4  = 0;

% high-wavenumber filter switch
% (set equal to 1 if you DON'T want filtering)
nofilter=0;

% planetary vorticity gradient
beta = 1;

% wind stress forcing
F = 1e-03;


% final time
tfin=30/mu;

% time-step
dt = 0.05;


% plot fields every tplot
tplot = 1/mu;
Nplot=round(tplot/dt);


% create x-y grid
dx=Lx/Nx;x=0:dx:Lx-dx;x=x.';
dy=Ly/Ny;y=0:dy:Ly-dy;y=y.';
[ X, Y]=meshgrid( x, y);

% create kx-ky wavenumber grid and operators
kx=2*pi/Lx*[0:Nx/2-1 -Nx/2:-1];kx=kx.';
ky=2*pi/Ly*[0:Ny/2-1 -Ny/2:-1];ky=ky.';
[KX,KY]=meshgrid(kx,ky);

ksq=KX.^2+KY.^2;
del2=-ksq;

% set (kx,ky)=1 component =1 so that del2
% can be used for invert laplacian
del2(KX==0&KY==0)=1;


% this will be \bar{zeta}(x,y)
zhatMean = zeros(Ny,Nx);

% the time after which \bar{zeta} will be calculated
tMean = 10/mu;
  NtMean = round(tMean/dt);
iterMean = 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HIGH-WAVENUMBER FILTER
% a is chosen so that the energy at the largest nondim
% wavenumber K*dx be zero whithin machine double precision
s=4;
Kmax = Ny/2; Kmax_s=Kmax*dy;
kcut = 2/3*Kmax;kcut_s=kcut*dy;
a = -log(1e-15)/(Kmax_s-kcut_s)^s * dy^s;
K=sqrt(KX.^2+KY.^2);

FILTER = 1*ones(Ny,Nx).*abs(K<=Ny/3) + exp(-a*(K-kcut).^s).*abs(K>Ny/3);
FILTER(KX==0&KY==0)=0;
if nofilter==1, FILTER = ones(Ny,Nx);end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% topographic potential vorticity eta=-f0*h/H
load topo_monoscale_v1;Nx_eta=256;Ny_eta=256;

[nx_eta,ny_eta]=size(eta);
if Nx<nx_eta || Ny<ny_eta
    error('to load this topography use higher resolution');
end

selx=[1:Nx_eta/2 Nx-Nx_eta/2+1:Nx];
sely=[1:Ny_eta/2 Ny-Ny_eta/2+1:Ny];

% if Nx>Nx_eta or Ny>Ny_eta then fit eta to the larger grid
etahat=zeros(size(KX));
etahat(sely,selx)=fft2(eta)*(Nx*Nx/(Nx_eta*Ny_eta));
eta = real(ifft2(etahat));

% <eta^2>
etarms = sqrt(mean(eta(:).^2));

% normalize eta so that
eta = eta/etarms;
etahat=fft2(eta);

test_norm_eta=sum(abs(etahat(:)).^2)/(Nx*Ny)^2;
if abs(test_norm_eta-1)>1e-10
    error('eta was not normalized correctly');
end

% calculate eta_x and eta_y
detadx = real(ifft2(1i*KX.*etahat));
detady = real(ifft2(1i*KY.*etahat));

% calculate r.m.s. of gradient of eta
gradetarms = sqrt(mean(detadx(:).^2+detady(:).^2));

% calculate ell_eta
ell = etarms/gradetarms;


% define the time vector
T = 0:dt:tfin;

% define vectors for saving variables
EU=0*T;         % E_U
Epsi=0*T;       % E_psi
Ut=0*T;         % U
Qpsi =0*T;      % Q_psi
formstress=0*T; % formstress


% initial flow fields: zhat and U
zhat = zeros(Ny,Nx);
   U = 0;

% initial values for formstress, U, E_U, E_psi and Q_psi
  EU(1) =  0.5*U^2;
Epsi(1) = -0.5*sum(abs(zhat(:)).^2./del2(:))/(Nx*Ny)^2;
Qpsi(1) =  0.5*sum(abs(zhat(:)).^2)/(Nx*Ny)^2;
  Ut(1) = U;
formstress(1) = real(mean(conj(zhat(:)./del2(:)).*(1i*KX(:).*etahat(:)))) / (Nx*Ny);


% calculate coefficients needed for the ETDRK4 time-step
[eL,eL2,Qfac,fu,fab,fc] = ETDRK4_coeffs(beta,mu,nu2,nu4,ksq,KX,del2,dt,Nx,Ny);


for it=2:length(T)


    % calculate \bar{zeta}
    if it>NtMean
      zhatMean = zhatMean*(1-1/iterMean) + zhat/iterMean ;
      iterMean = iterMean + 1;
    end

    % time-step U,zhat
    [zhatnew,Unew] = time_step_ETDRK4(zhat,U,eL,eL2,Qfac,fu,fab,fc,dt);
    U = Unew;
    zhat = zhatnew.*FILTER;

    % save U(t)
    Ut(it)=U;

    % save E_U(t)
    EU(it) = 0.5*U^2;

    % save E_psi(t)
    Epsi(it) = -.5*sum(abs(zhat(:)).^2./del2(:))/(Nx*Ny)^2;

    % save formstress(t)
    formstress(it) = real(mean(conj(zhat(:)./del2(:)).*(1i*KX(:).*etahat(:)))) / (Nx*Ny);

    % save Q_psi(t)
    Qpsi(it) = 0.5*sum(abs(zhat(:)).^2)/(Nx*Ny)^2;

    % every Nplot time-steps plot fields
    if rem(it,Nplot)==1||it==length(T)

        % calculate fields in physical space
        psihat = zhat./del2;
        zeta=real(ifft2(zhat));
        psi = real(ifft2(psihat));
        zetaMean = real(ifft2(zhatMean));
        u = real(ifft2(-1i*KY.*psihat));
        v = real(ifft2(+1i*KX.*psihat));

        % check the CFL
        cfl = sqrt(max((U+u(:)).^2+v(:).^2))*dt/dx;
        display(['t=' num2str(T(it)) '  cfl=' num2str(cfl,'%1.3f')]);
        if cfl>.8
            error('cfl > 0.8, reduce dt');
        end

        % plot fields
        plot_figs
    end

end
