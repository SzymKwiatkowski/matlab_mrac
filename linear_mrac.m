close all; clear all;

%% Drone Parameter Values
global Ka Km m Ix Iy Iz g l;
Ct = 0.0107;
Cq = Ct*sqrt(Ct/2);
Rr = 33/1000;    % rotor radius
RA = pi*Rr^2;     % rotor radius
rho = 1.184;    % density of air
Ka = Ct*rho*RA*Rr^2;
Km = Cq*rho*RA*Rr^3;
g = 9.81;
Ix = 0.0686e-3;
Iy = 0.092e-3;
Iz = 0.1366e-3;
l = 0.0624;        %Distance from rotor to the center of Drone
m = 0.068;

 %% Linearized System arround the equilibrium point (Hover)
global K A B;

O6 = zeros(6);
I6 = eye(6);
Psi = 0;    % Phi angule chosen as Eq point at the hover.

Phi = [0, 0, 0, -g*sin(Psi), -g*cos(Psi), 0;
       0, 0, 0, g*cos(Psi), -g*sin(Psi), 0;
       0, 0, 0, 0,           0,          0;
       0, 0, 0, 0,           0,          0;
       0, 0, 0, 0,           0,          0;
       0, 0, 0, 0,           0,          0];
   
A = [O6,  I6;       % Jacobian Matrix A
     Phi, O6];
 
Ac = [Phi, O6];     % Part of A that is Controllable

O84 = zeros(8,4);

Delta = [Ka/m,      Ka/m,           Ka/m,       Ka/m;
         0,         -Ka*l/Ix,       0,          Ka*l/Ix;
         Ka*l/Iy,   0,              -Ka*l/Iy,   0;
         Km/Iz,     -Km/Iz,         Km/Iz,      -Km/Iz];
       
B = [O84;       % Jacobian Matrix B
     Delta];
 
%% Model Mismatch (the actual plant)
global A_act B_act
A_act = A;
m_act = 1.5 * m;
l_act = 1.5 * l;
Delta_act = [   Ka/m_act,     Ka/m_act,     Ka/m_act,    Ka/m_act;
               0, -Ka*l_act/Ix,        0, Ka*l_act/Ix;
         Ka*l_act/Iy,        0, -Ka*l_act/Iy,       0;
           Km/Iz,   -Km/Iz,    Km/Iz,  -Km/Iz];

B_act = [O84;
         Delta_act];
 
%% Pole Placement
global Kr_lin_ctr;

desired_poles = -linspace(1,12,12);
K = place(A,B,desired_poles);
Kr_lin_ctr = B\(A-B*K);

% Our calculated Gain K to stabilize the system.
array2table(K)

% Eigenvalues of the closed loop system.
EigenValues = array2table(eig(A - B*K))

%% Set the desired State Values (Reference Signal that the 12 states should track)
global r;
r = [1,1,1,0,0,0,0,0,0,0,0,0]';     % the desired state values


%% Linear Adaptive Controller & Linearized Plant (Mismatched Parameters)
tspan = [0,5];
x0 = zeros(88,1);
% initial conditions: x(0) = xm(0) = 0
% use the ideal linearized plant to obtain kx(0) and kr(0)
temp_K = -K';
x0(25:72) = [temp_K(1,:) temp_K(2,:) temp_K(3,:) temp_K(4,:) temp_K(5,:) temp_K(6,:) temp_K(7,:) temp_K(8,:) temp_K(9,:) temp_K(10,:) temp_K(11,:) temp_K(12,:)]';
x0(73) = 1;
x0(78) = 1;
x0(83) = 1;
x0(88) = 1;
[t, x] = ode45(@AC_LinearMismatchedModel, tspan, x0);
x_adap_act = x(:,1:12);
xm_adap_act = x(:,13:24);
%% Plots
for i = 1:12
    figure;
    hold on
    plot(t,x_adap_act(:,i));
    plot(t,xm_adap_act(:,i));
    legend(['x-adaptive_{' num2str(i) '}(t)'], ['xm_{' num2str(i) '}(t)']);
    xlabel('t');
    title('Adaptive & Linear Controller with Imperfect Plant: x-adaptive(t), x-linear(t), xm(t)');
    grid on;
    hold off
end
