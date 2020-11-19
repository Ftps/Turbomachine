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
t_a1 = (psi/2 + 1 - r)/phi;
t_b1 = t_a1 - 1/phi;
t_b2 = (phi*t_a1 - 1 - psi)/phi;
t_a2 = t_b2 + 1/phi;

a1 = atand(t_a1);
b1 = atand(t_a2);
a2 = atand(t_a2);
b2 = atand(t_b2);

% Meridional rotational speed and inlet velocity
U = sqrt(w/psi);
vc = phi*U;
v1 = vc/cosd(a1);

% Mean, top and bottom radii (inlet), and angular speed of the rotor
T1 = T01 - (v1^2)/(2*cp);
Umax = 0.9*sqrt(k*R*T1);	% Limiting condition for the top blade height
rho1 = (P01/(R*T01))*(T1/T01)^(1/(k-1));

options = optimoptions('fsolve','Display','iter');
f = @(x) [U - x(3)*(x(1)+x(2))/2, Umax - x(3)*x(2), rho1*pi*(x(2)^2-x(1)^2)*vc - m];
x = fsolve(f, [0.5, 1, 356], options);

r1 = x(1);
r2 = x(2);
ww = x(3);
rm = 0.5*(r1+r2);
