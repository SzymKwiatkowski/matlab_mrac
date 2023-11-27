%% Linear Controller & Ideal Linearized Plant
function xdot = LC_LinearModel(t,X)
    global K A B Kr_lin_ctr r;
    ut = -K*X - Kr_lin_ctr*r;
    xdot = A*X + B*ut;
end
