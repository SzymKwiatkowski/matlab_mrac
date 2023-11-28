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

%% Check if controllable
Co = ctrb(A,B);
unco = length(A) - rank(Co)

%% Model Mismatch (the actual plant)
global A_act B_act
A_act = A;
m_act = 1 * m;
l_act = 1 * l;
Delta_act = [   Ka/m_act,     Ka/m_act,     Ka/m_act,    Ka/m_act;
               0, -Ka*l_act/Ix,        0, Ka*l_act/Ix;
         Ka*l_act/Iy,        0, -Ka*l_act/Iy,       0;
           Km/Iz,   -Km/Iz,    Km/Iz,  -Km/Iz];

B_act = [O84;
         Delta_act];

%% Check if controllable
Co = ctrb(A_act,B_act);
unco = length(A) - rank(Co)

%% Pole Placement
global Kr_lin_ctr;

desired_poles = -linspace(0.25,2.0,12);
K = place(A,B,desired_poles);
Kr_lin_ctr = B\(A-B*K);

% Our calculated Gain K to stabilize the system.
array2table(K)

% Eigenvalues of the closed loop system.
EigenValues = array2table(eig(A - B*K));

%% Set the desired State Values (Reference Signal that the 12 states should track)
global r;
r = [1,1,1,0,0,0,0,0,0,0,0,0]';     % the desired state values

%% Linear Controller & Nonlinear Plant
tspan = [0,3.5];
x0 = zeros(1,12);
[t_lc_nlp, x_lc_nlp] = ode45(@LC_NonlinearModel, tspan, x0);

%% Plots
figure;
subplot(4,1,1);
hold on
plot(t_lc_nlp, x_lc_nlp(:,1));
plot(t_lc_nlp, x_lc_nlp(:,2));
plot(t_lc_nlp, x_lc_nlp(:,3));
title('State Trajectories of Fixed-gain Linear Controller with Nonlinear Plant');
legend('x_{1}(t)','x_{2}(t)','x_{3}(t)');
xlabel('t');
ylabel('Magnitude');
grid on

subplot(4,1,2);
hold on
plot(t_lc_nlp, x_lc_nlp(:,4));
plot(t_lc_nlp, x_lc_nlp(:,5));
plot(t_lc_nlp, x_lc_nlp(:,6));
legend('x_{4}(t)','x_{5}(t)','x_{6}(t)');
xlabel('t');
ylabel('Magnitude');
grid on

subplot(4,1,3);
hold on
plot(t_lc_nlp, x_lc_nlp(:,7));
plot(t_lc_nlp, x_lc_nlp(:,8));
plot(t_lc_nlp, x_lc_nlp(:,9));
legend('x_{7}(t)','x_{8}(t)','x_{9}(t)');
xlabel('t');
ylabel('Magnitude');
grid on

subplot(4,1,4);
hold on
plot(t_lc_nlp, x_lc_nlp(:,10));
plot(t_lc_nlp, x_lc_nlp(:,11));
plot(t_lc_nlp, x_lc_nlp(:,12));
legend('x_{10}(t)','x_{11}(t)','x_{12}(t)');
xlabel('t');
ylabel('Magnitude');
grid on

hold off


