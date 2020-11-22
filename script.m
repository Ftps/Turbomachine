% Air/Physical properties
bar = 1e5;
k = 1.4;
cp = 1000;
R = 8.3144/(28.9647e-3);

% Given values
m_cor = 54.4;
pp_total = 23;
n_ad = 0.86;
n_pl = 0.906;
T01 = 304.8;
P01 = 0.60496*bar;

% Real mass flow rate, stage number and stage pressure ratio
m = m_cor*(P01/101325)/sqrt(T01/288.15);
n_stages = ceil(9*log(pp_total)/log(14));
pp = pp_total^(1/n_stages);

% Stage specific work and isentropic efficiency
w = cp*T01*(pp^((k-1)/(n_pl*k)) - 1);
n_is = (pp^((k-1)/k) - 1)/(pp^((k-1)/(n_pl*k)) - 1);

% Related values from the Smith charts
r = 0.5;
psi = 0.4;
phi = 0.9;

% Velocity triangle angles
t_a1 = (-psi/2-r+1)/phi;
t_b1 = t_a1 - 1/phi;
t_b2 = (phi*t_a1 - 1 + psi)/phi;
t_a2 = t_b2 + 1/phi;

a1 = atand(t_a1);
b1 = atand(t_b1);
a2 = atand(t_a2);
b2 = atand(t_b2);

% Meridional rotational speed and velocities
U = sqrt(w/psi);
vc = phi*U;
v1 = vc/cosd(a1);
v2 = vc/cosd(a2);

% Mean, top and bottom radii (inlet), and angular speed of the rotor
T1 = T01 - (v1^2)/(2*cp);
Umax = 0.9*sqrt(k*R*T1);	% Limiting condition for the top blade height
rho1 = (P01/(R*T01))*(T1/T01)^(1/(k-1));

% Constraints for the radii and rotacional speed:
% 1) At the tip, speed is equal to Umax:
%		Umax = r2*ww;
% 2) At the mid point, velocity is equal to U:
%		U = 0.5*(r1+r2)*ww;
% 3) Mass flow rate must be the calculated value:
%		rho1*pi*(r2^2 - r1^2)*vc = m
%
% Transform this into a problem f(x) = 0 and solve
%
%					/ U - ww*0.5*(r1+r2)
%	f(r1, r2, ww) = | Umax - ww*r2
%					\ rho1*pi*(r2^2-r1^2)*vc - m

options = optimoptions('fsolve','Display','iter');
f = @(x) [U - x(3)*(x(1)+x(2))/2, Umax - x(3)*x(2), rho1*pi*(x(2)^2-x(1)^2)*vc - m];
x = fsolve(f, [0.5, 1, 356], options);

r1 = x(1);
r2 = x(2);
ww = x(3);
rm = 0.5*(r1+r2);

% Allocate space for the properties at all stage positions
p = zeros(2*n_stages+1, 1);
p0 = zeros(2*n_stages+1, 1);
T = zeros(2*n_stages+1, 1);
T0 = zeros(2*n_stages+1, 1);
R1 = zeros(2*n_stages+1, 1);
R2 = zeros(2*n_stages+1, 1);

% Vector for plotting
xx = (1:(2*(n_stages+1)))./2;
xx = xx(2:end);

% Cycle that moves through all stage positions
for i = 1:(2*n_stages+1)
	% Stage one inlet
	if i == 1
		p0(1) = P01;
		T0(1) = T01;
		R1(1) = r1;
		R2(1) = r2;
		T(1) = T1;
		p(1) = P01*(T1/T01)^(k/(k-1));
	% Position between rotor and stator
	elseif mod(i,2) == 0
		p0(i) = pp*p0(i-1);
		T0(i) = T0(i-1)*(pp^((k-1)/(k*n_pl)));
		T(i) = T0(i) - (v2^2)/(2*cp);
		p(i) = p0(i)*(T(i)/T0(i))^(k/(k-1));
		
		rho1 = p(i)/(R*T(i));
		f = @(x) [rm - 0.5*(x(2) + x(1)), rho1*pi*(x(2)^2-x(1)^2)*vc - m];
		x = fsolve(f, [R1(i-1), R2(i-1)], options);
		R1(i) = x(1);
		R2(i) = x(2);
	% Position between stages
	else
		p0(i) = p0(i-1);
		T0(i) = T0(i-1);
		T(i) = T0(i) - (v1^2)/(2*cp);
		p(i) = p0(i)*(T(i)/T0(i))^(k/(k-1));
		
		rho1 = p(i)/(R*T(i));
		f = @(x) [rm - 0.5*(x(2) + x(1)), rho1*pi*(x(2)^2-x(1)^2)*vc - m];
		x = fsolve(f, [R1(i-1), R2(i-1)], options);
		R1(i) = x(1);
		R2(i) = x(2);
	end
end
