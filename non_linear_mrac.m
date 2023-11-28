close all; clear all; clc;

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

%% Nonlinear Adaptive Controller & Nonlinear Plant
global lambda gamma_x gamma_r gamma_alpha;
lambda = 0.9; % the parameter in the coord transform

% Some Preparation
global A_ref B_ref K_ref Kr_nonlin_ctr r_pos Am Bm P;
A_ref = [zeros(3,3) , eye(3) ; zeros(3,3) , zeros(3,3)];
B_ref = zeros(6,4);
B_ref(4,1) = -lambda*Ka*l/Iy;
B_ref(4,3) = lambda*Ka*l/Iy;
B_ref(5,2) = -lambda*Ka*l/Ix;
B_ref(5,4) = lambda*Ka*l/Ix;
B_ref(6,:) = [-Ka/m , -Ka/m, -Ka/m, -Ka/m];
desired_poles = [-1,-2,-3,-4,-4,-3];
K_ref = place(A_ref,B_ref,desired_poles);
Kr_nonlin_ctr = B_ref\(A_ref-B_ref*K_ref);

% the reference model
Am = A_ref - B_ref*K_ref;
Bm = B_ref;
% matrix Q , P
Q = 180*eye(6);
P = lyap(Am',Q);
% adaptation rates
gamma_x = 80*eye(6);
gamma_r = 80*eye(4);
gamma_alpha = 80*eye(6);

% The desired values for the 6 position states
% x1,x2,x3,x7,x8,x9 (in the original coords)
r_pos = [0;0;1;0;0;0];

% start the simulation
tspan = [0, 0.8];
x0 = zeros(88,1);
% initial conditions: x(0) = xm(0) = 0
% convert xm(0) to the new transformed coordinates
x0(15) = -lambda;
x0(85) = -lambda;
% use the position state dynamics in the new coords to obtain kx(0) and kr(0)
temp_K = -K_ref';
x0(19:42) = [temp_K(1,:) temp_K(2,:) temp_K(3,:) temp_K(4,:) temp_K(5,:) temp_K(6,:)]';
x0(43) = 1;
x0(48) = 1;
x0(53) = 1;
x0(58) = 1;

% calculate alpha_hat(0)
temp1 = [-Ka/m , -Ka/m , -Ka/m, -Ka/m;
         0 , -lambda*Ka*l/Ix , 0 , lambda*Ka*l/Ix;
         -lambda*Ka*l/Iy , 0 , lambda*Ka*l/Iy , 0];
temp2 = [lambda , g , 0 , 0 , 0 , 0;
         0 , 0 , lambda*(-Ix+Iy-Iz)/Ix , g , 0 , 0;
         0 , 0 , 0 , 0 , -lambda*(-Ix+Iy+Iz)/Iy , -g];
alpha_star = (temp1\temp2)';
x0(59:82) = [alpha_star(1,:) , alpha_star(2,:) , alpha_star(3,:) , alpha_star(4,:) , alpha_star(5,:) , alpha_star(6,:)]';
%% Calculate result
[t, x] = ode45(@AC_NonlinearModel, tspan, x0);

%% Plots
xt = x(:,1:12);
xm = x(:,13:18);
x_prime = x(:,83:88);

for i = 1:6
    figure;
    hold on
    plot(t,x_prime(:,i));
    plot(t,xm(:,i));
    xlabel('t');
    legend(['x_{' num2str(i) '}(t)'],['xm_{' num2str(i) '}(t)']);
    title('x(t) and xm(t) for the 6 position states (in the transformed coordinates)');
    grid on
    hold off
end

for i = 1:12
    figure;
    plot(t,xt(:,i));
    xlabel('t');
    legend(['x_{' num2str(i) '}(t)']);
    title('x(t) in the original coordinates');
    grid on;
end
