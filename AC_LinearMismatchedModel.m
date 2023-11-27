%% Linear Adaptive Controller(MRAC MIMO) & Linearized Plant (mismatched parameters)
function xdot = AC_LinearMismatchedModel(t,x)
global K Kr_lin_ctr A B A_act B_act r;       % x: [x xm kx(1-4,:) km(1-4,:)]

disp(t);
Am = A - B*K;
Bm = B;
Q = 600*eye(12);          % 600
P = lyap(Am',Q);
gamma_x = 0.005*eye(12);    % 0.005
gamma_r = 0.005*eye(4);     % 0.005
rt = -Kr_lin_ctr * r;

x = reshape(x,[1,88]);
xt = x(1:12)';
xm = x(13:24)';
Kx = [x(25:28);x(29:32);x(33:36);x(37:40);x(41:44);x(45:48);x(49:52);x(53:56);x(57:60);x(61:64);x(65:68);x(69:72)];
Kr = [x(73:76) ; x(77:80) ; x(81:84) ; x(85:88)];
ut = Kx' * xt + Kr' * rt;

e = xt - xm;

Kx_dot = -gamma_x * xt * e' * P * B_act;
Kr_dot = -gamma_r * rt * e' * P * B_act;

xdot = zeros(88,1);
xdot(1:12) = A_act*xt + B_act*ut;
xdot(13:24) = Am*xm + Bm*rt;
xdot(25:72) = [Kx_dot(1,:) Kx_dot(2,:) Kx_dot(3,:) Kx_dot(4,:) Kx_dot(5,:) Kx_dot(6,:) Kx_dot(7,:) Kx_dot(8,:) Kx_dot(9,:) Kx_dot(10,:) Kx_dot(11,:) Kx_dot(12,:)]';
xdot(73:88) = [Kr_dot(1,:) Kr_dot(2,:) Kr_dot(3,:) Kr_dot(4,:)]';
end
