% problem setup script
% higher perm

oneyear = 60*60*24*365.25; % (s) 1 yr

% parameters are stored in data structure M

% along-strike discretization and fluid transport
M.beta = 1e-8; % pore+fluid compressiblity (Pa^-1)
M.phi = 1e-1; % porosity
M.k = 5e-13; % permeability (m^2)
M.mu = 1e-3; % fluid viscosity (Pa s)
M.c = M.k/(M.beta*M.phi*M.mu); % hydraulic diffusivity (m^2/s)
M.L = 60e3; % along-strike distance (m)
refine = 2; % mesh refinement factor (refine=1 gives 1 km grid spacing)
M.nx = 60*refine; % number of intervals, nx+1 grid points including endpoints
M.h = M.L/M.nx; % grid spacing
M.x = [0:M.nx]'*M.h-M.L/2; % nx+1 grid points (column vector)
M.w = 15; % fault zone width (m)
M.H = (2.7-1.6)*1e3; % down-dip distance (m)
M.i = M.nx/2+1+[-8 -1 8]*refine; % indices of injector(s)
M.xi = M.x(M.i); % location of injector(s)
M.Q0 = [1/2.5 0.5/2.5 3/4]*1e6/oneyear; % volumetric injection rate (m^3/s)
M.ti = [1.5 1.5 0]*oneyear; % start time of injectors (s)
M.tstop = [Inf Inf Inf]*oneyear; % stop time of injectors (s)
M.q0 = M.Q0/(M.H*M.w); % linear injection rate (m/s)

% rate-and-state parameters
M.f0 = 0.6; % reference friction coefficient
M.V0 = 1e-6; % reference slip velocity
M.a = 0.01; % direct effect parameter
M.b = 0.009; % state evolution parameter
M.dc = 100e-6; % (m) state evolution distance

% elasticity parameters
rho = 2.54; % density (g/cm^3)
c = 2.5; % S-wave speed (km/s)
M.eta = rho*c/2; % radiation-damping coefficient (MPa*s/m)

% initial effective normal stress
M.N = 15; % (MPa)

% stiffness (elastic stress change per unit slip)
M.K = 16/(1-0.26)/(M.H*1e-3); % (MPa/m)

% loading:
% shear stress increases at constant rate in absence of slip
% pressure increases at constant rate
M.dtaudt = 0; % (MPa/s)
M.dpdt = 0; % pressurization rate (MPa/s)

% initial conditions
D0 = zeros(M.nx+1,1); % slip (m)
M.tau0 = 0.545*M.N; % shear stress (MPa)
Psi0 = 0.7*ones(M.nx+1,1); % state variable
p0 = zeros(M.nx+1,1); % pressure (MPa)

% time
tmax = 3*oneyear;

% adaptive time-stepping options

dispAcceptReject = false;
plotSolDuringSim = false;

tol = 1e-4; % tolerance
dt = 1; % initial time step (s)
dtmax = 1e6; % maximum time step (s)
safety = 0.9; % safety factor
