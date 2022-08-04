%% Generates Table 1

close all; clear; clc

configuration = 'outward'; d = 1;
% configuration = 'outward'; d = 2;
% configuration = 'outward'; d = 3;
% configuration = 'inward'; d = 1;
% configuration = 'inward'; d = 2;
% configuration = 'inward'; d = 3;

%% Material and laser flash parameters
Qinf = 7000; % total heat absorbed
l0 = 0.001; % inner surface at x = l0
l1 = 0.003; % outer surface at x = l1
k = 222; % thermal conductivity
rho = 2700; % density
cp = 896; % specific heat capacity
N = 1000; % number of temperature rise values (excluding t = 0)
tN = 0.1; % end time
beta = 0.001; % exponential pulse parameter (peak occurs at t = beta)
q = @(t) Qinf*t.*exp(-t/beta)/beta^2; % exponential pulse
dt = tN/N; % time step duration
t = (0:dt:tN)'; % discrete times
alpha = k/(rho*cp); % target value of thermal diffusivity

%% Finite volume method parameters
Nr = 501; % number of nodes
h = (l1-l0)/(Nr-1); % node spacing
x = linspace(l0,l1,Nr); % location of nodes
xw(1) = x(1); xw(2:Nr) = (x(1:Nr-1)+x(2:Nr))/2; % west boundaries
xe(1:Nr-1) = xw(2:Nr); xe(Nr) = x(Nr); % east boundaries
T0 = zeros(1,Nr); % initial temperature rise
AbsTol = 1e-12; % absolute error tolerence
RelTol = 1e-12; % relative error tolerence

% Sparsity pattern of Jacobian
e = ones(Nr,1); JPat = spdiags([e e e],-1:1,Nr,Nr);
options = odeset('RelTol',RelTol,'AbsTol',AbsTol,'JPattern',JPat);

%% Computes entries in Table 1

% Solve heat transfer model
[~,T] = ode15s(@(t,T) Gfunc(t,T,d,alpha,l0,l1,h,x,xw,xe,rho,cp,Nr,configuration,q),t,T0,options);
if isequal(configuration,'outward')
    Tdata = T(:,end); % temperature rise at outer surface
    Tinf = d*l0^(d-1)*Qinf/(rho*cp*(l1^d-l0^d)); % steady state temperature rise
elseif isequal(configuration,'inward')
    Tdata = T(:,1); % temperature rise at inner surface
    Tinf = d*l1^(d-1)*Qinf/(rho*cp*(l1^d-l0^d)); % steady state temperature rise
end

% Estimate of thermal diffusivity
alpha_tilde = thermal_diffusivity(d,t,Tdata,Tinf,l0,l1,beta);
display(alpha_tilde);

% Signed relative percent error
epsilon = (alpha - alpha_tilde)/alpha * 100;
display(epsilon)
