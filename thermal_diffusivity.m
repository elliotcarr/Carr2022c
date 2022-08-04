function alpha = thermal_diffusivity(d,t,Tr,Tinf,l0,l1,beta)
% Calculate thermal diffusivity according to equations (18) and (24) with
% trapezoidal rule approximations (25) and (26)

% Inputs
%   d: dimension (d = 1,2 or 3)
%   t: discrete times (t0,...,tN)
%   Tr: temperature rise values at discrete times (Tr0,...,TrN)
%   Tinf: steady state temperature rise
%   l0: inner surface location 
%   l1: outer surface location (l1 > l0)
%   beta: exponential pulse parameter (peak occurs at t = beta)
% Outputs
%   alpha: thermal diffusivity

IT = trapzrule(t,1-Tr/Tinf);
Iq = 2*beta;
if d == 1
I = l1 - l0;
elseif d == 2
    I = log(l1/l0);
elseif d == 3
    I = 1/l0 - 1/l1;
end
alpha = (l1^(d+2)-(d+2)*l0^d*l1^d*I-l0^(d+2)) / (2*(d+2)*(l1^d-l0^d)*(IT - Iq));

% syms s;
% alpha = (l1^(d+2)-(d+2)*l0^d*l1^d*double(int(s^(1-d),l0,l1))-l0^(d+2)) ...
%     / (2*(d+2)*(l1^d-l0^d)*(IT - Iq));
end

function integral = trapzrule(t,f)
% Trapezoidal rule

integral = sum((t(2:end)-t(1:end-1)).*((f(1:end-1)+f(2:end))/2));

end