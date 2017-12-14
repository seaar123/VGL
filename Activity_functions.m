function [outs,D_outs]=Activity_functions(ins,knd)
if knd == 1
    %knd=1 linear ( identity )
    outs=ins;% y
    D_outs=1;% dy/dt
elseif knd == 2
    %knd=2 logistic sigmmid
    outs=1/(1+exp(-ins));% y
    D_outs=outs*(1-outs);% dy/dt
elseif knd == 3
    %knd=3 hyperbolic tangent ( tanh() )
    outs=(exp(ins)-exp(-ins))/(exp(ins)+exp(-ins));% y
    D_outs=1-(tanh(ins))^2;% dy/dt
elseif knd == 4
    %knd=4 gaussion
    outs=(1-exp(-ins))/(1+exp(-ins));% y
    D_outs=0.5*(1-(tanh(ins))^2);% dy/dt
end