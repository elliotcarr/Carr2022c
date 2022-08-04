function G = Gfunc(t,T,d,alpha,l0,l1,h,x,xw,xe,rho,cp,Nx,configuration,q)
% Finite volume discretisation in space is expressible in the form:
%           dT/dt = G(t,T)
% where T = (T1,...,TNx) contains the temperature rise values at nodes i =
% 1,...,Nx. Gfunc evalates the right-hand side function G(t,T).

% Inputs
%   t: time
%   T: temperature rise nodal values (T1,...,TNx)
%   d: dimension (d = 1,2 or 3)
%   alpha: thermal diffusivity
%   l0: inner surface location
%   l1: outer surface location (l1 > l0)
%   h: uniform node spacing
%   x: node locations
%   xw: west control volume boundaries
%   xe: east control volume boundaries
%   rho: density
%   cp: specific heat capacity
%   Nx: number of nodes
%   configuration: 'outward' or 'inward'
%   q: heat flux function
% Outputs
%   G: vector-valued function G evaluated at t and T

G = zeros(Nx,1);

% Boundary nodes
if isequal(configuration,'outward') % heat pulse applied at inner surface (x = l0)
    G(1) = 2*alpha*xe(1)^(d-1)*(T(2)-T(1))/(x(1)^(d-1)*h^2) + 2*l0^(d-1)*q(t)/(x(1)^(d-1)*h*rho*cp);
    G(Nx) = 2*alpha*xw(Nx)^(d-1)*(T(Nx-1)-T(Nx))/(x(Nx)^(d-1)*h^2);
elseif isequal(configuration,'inward') % heat pulse applied at outer surface (x = l1)
    G(1) = 2*alpha*xe(1)^(d-1)*(T(2)-T(1))/(x(1)^(d-1)*h^2);
    G(Nx) = 2*l1^(d-1)*q(t)/(x(Nx)^(d-1)*h*rho*cp) - 2*alpha*xw(Nx)^(d-1)*(T(Nx)-T(Nx-1))/(x(Nx)^(d-1)*h^2);
end

% Interior nodes
for i = 2:Nx-1
    G(i) = alpha*(xe(i)^(d-1)*T(i+1)-(xe(i)^(d-1)+xw(i)^(d-1))*T(i)+xw(i)^(d-1)*T(i-1))/(x(i)^(d-1)*h^2);
end

end