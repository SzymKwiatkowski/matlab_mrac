%% Linear Controller & Linearized Plant (mismatched parameters)
function xdot = LC_LinearMismatchedModel(t,x)
    global K A_act B_act Kr_lin_ctr r;
    ut = -K*x - Kr_lin_ctr*r;
    xdot = A_act*x + B_act*ut;
end