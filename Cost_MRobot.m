function [U,dU_dx,dU_du]=Cost_MRobot(x,u)
tau=u';
x=x';
%--------------------------------------------------------------------------
% Cost
Q=[1 0;0 1];
R=[1 0;0 1];
U=(x)'*Q*(x)+tau'*R*tau;
dU_dx=2*Q*(x);
dU_du=2*R*tau;
