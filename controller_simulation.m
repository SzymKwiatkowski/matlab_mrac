%% Adaptive Project - Quadcopter

%% Drone Parameter Values
%close all; clear all;
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

%% Linear Controller & Ideal Linearizd Plant
tspan = [0, 5];
X0 = [0;0;0;0;0;0;0;0;0;0;0;0]; % [x, y, z, phi, theta, psi, xdot, ydot, zdot, p, q, r] initial conditions.

[t, X] = ode45(@LC_LinearModel, tspan, X0); % Here X is actually Xtilde, where Xtilde = X - Xeq. 

%Graph of the states of Xtilde, where Xtilde = X - Xeq..
titlelegend = ["Positions in x, y, z", "Euler angules: \phi, \theta, \psi.", "Velocities in x, y, z.", "Body angular rates: p, q, r."];
i = 1;
for j = 1:4
    figure; hold on;
    plot(t,X(:,i)); plot(t,X(:,i+1)); plot(t,X(:,i+2));
    legend(['x_{' num2str(i) '}(t)'],['x_{' num2str(i+1) '}(t)'],['x_{' num2str(i+2) '}(t)']);
    xlabel('t');
    ylabel('Magnitude of the state $\tilde{X}$', 'Interpreter','latex')
    title(titlelegend(j));
    grid on;
    i = i + 3;
end

Xeq = [0,0,3,0,0,Psi,0,0,0,0,0,0]; % [x, y, z, phi, theta, psi, xdot, ydot, zdot, p, q, r] at eq point.

%Graph of the states of X, where X = Xtilde + Xeq.
Xreal = zeros(length(X(:,1)),length(X(1,:)));
for i = 1:length(X(:,1))
    Xreal(i,:) = X(i,:) + Xeq;
end

i = 1;
for j = 1:4
    figure; hold on;
    plot(t,Xreal(:,i)); plot(t,Xreal(:,i+1)); plot(t,Xreal(:,i+2));
    legend(['x_{' num2str(i) '}(t)'],['x_{' num2str(i+1) '}(t)'],['x_{' num2str(i+2) '}(t)']);
    xlabel('t');
    ylabel('Magnitud of the actual state X', 'Interpreter','latex')
    title(titlelegend(j));
    grid on;
    i = i + 3;
end

%% Linear Controller & Linearized Plant (Mismatched Parameters)
tspan = [0,5];
x0 = zeros(1,12);
[t_lin_act, x_lin_act] = ode45(@LC_LinearMismatchedModel, tspan, x0);

% Graph of the states of Xtilde, where Xtilde = X - Xeq..
titlelegend = ["Positions in x, y, z", "Euler angules: \phi, \theta, \psi.", "Velocities in x, y, z.", "Body angular rates: p, q, r."];
i = 1;
for j = 1:4
    figure; hold on;
    plot(t,X(:,i)); plot(t,X(:,i+1)); plot(t,X(:,i+2));
    legend(['x_{' num2str(i) '}(t)'],['x_{' num2str(i+1) '}(t)'],['x_{' num2str(i+2) '}(t)']);
    xlabel('t');
    ylabel('Magnitude of the state $\tilde{X}$', 'Interpreter','latex')
    title(titlelegend(j));
    grid on;
    i = i + 3;
end

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
    plot(t_lin_act,x_lin_act(:,i));
    legend(['x-adaptive_{' num2str(i) '}(t)'],['xm_{' num2str(i) '}(t)'],['x-linear_{' num2str(i) '}(t)']);
    xlabel('t');
    title('Adaptive & Linear Controller with Imperfect Plant: x-adaptive(t), x-linear(t), xm(t)');
    grid on;
    hold off
end

%% Linear Controller & Nonlinear Plant
tspan = [0,5];
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


%% Nonlinear Adaptive Controller & Nonlinear Plant
global lambda gamma_x gamma_r gamma_alpha;
lambda = 0.01; % the parameter in the coord transform

% Some Preparation
global A_ref B_ref K_ref Kr_nonlin_ctr r_pos Am Bm P;
A_ref = [zeros(3,3) , eye(3) ; zeros(3,3) , zeros(3,3)];
B_ref = zeros(6,4);
B_ref(4,1) = -lambda*Ka*l/Iy;
B_ref(4,3) = lambda*Ka*l/Iy;
B_ref(5,2) = -lambda*Ka*l/Ix;
B_ref(5,4) = lambda*Ka*l/Ix;
B_ref(6,:) = [-Ka/m , -Ka/m, -Ka/m, -Ka/m];
desired_poles = [-1,-2,-3,-4,-5,-6];       %linspace(5,10,6);
K_ref = place(A_ref,B_ref,desired_poles);
Kr_nonlin_ctr = B_ref\(A_ref-B_ref*K_ref);

% the reference model
Am = A_ref - B_ref*K_ref;
Bm = B_ref;
% matrix Q , P
Q = 300*eye(6);
P = lyap(Am',Q);
% adaptation rates
gamma_x = 50*eye(6);
gamma_r = 50*eye(4);
gamma_alpha = 50*eye(6);

% The desired values for the 6 position states
% x1,x2,x3,x7,x8,x9 (in the original coords)
r_pos = [0;0;1;0;0;0];

% start the simulation
tspan = [0,1.1];
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
%%
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




