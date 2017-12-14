function [tdot,MR_states]=Intgration_fun(t,h,MR_states)
FR=0; % =1 with friction =0 without friction
z=MR_states(1:5);%  x,y,theta,phiR, and phiL
y=MR_states(6:7);%  v ,and w
%--------------------------------------------------------------------------
%Robot parameters
mT=10;  %Total mass
mW=2;   %Mass of 1 wheel
r=0.05; %Wheel radius
b=0.4;  %half the robot width
d=0.1;  %CG offset from rear axle
Iyy=1;  %Wheel moment of Iertia
IT=5;   %Platform total moment of inertia
Fv=.5 *FR;  %Coefficient of Fiscous friction ( *0 means ignoring friction )
Fd=.8 *FR;  %Coefficient of Colomb friction
%--------------------------------------------------------------------------
%Robot States
xc=z(1);
yc=z(2);
th=z(3);
phR=z(4);
phL=z(5);
v=y(1);
w=y(2);
V=[v;w];
%--------------------------------------------------------------------------
%Robot Kinematics
S=[ cos(th)  -d*sin(th);
    sin(th)  d*cos(th);
    0        1;
    1/r      b/r;
    1/r     -b/r];
%--------------------------------------------------------------------------

xdot(1:5)=S*V;

phRdot=xdot(4);
phLdot=xdot(5);
%--------------------------------------------------------------------------

xdot(6:7)=y;

% z=[x;y;th;phiR; phiL]; z(1,:)=x  z(2,:)=y  z(3,:)=th  z(4,:)=phiR z(5,:)=phiL

g=@(t,z) xdot(1:5);
% y=[v;w]; y(1,:)=v  y(2,:)=w

f=@(t,y) xdot(6:7);

%--------------------------------------------------------------------------
%RK
%Updating time
tdot=t+h;
%Updating of y and z
ky1=f(t    ,y           );
kz1=g(t    ,z           );
ky2=f(t+(h/2),y+(h/2)*ky1);
kz2=g(t+(h/2),z+(h/2)*kz1);
ky3=f(t+(h/2),y+(h/2)*ky2);
kz3=g(t+(h/2),z+(h/2)*kz2);
ky4=f(t+h    ,y+ h   *ky3);
kz4=g(t+h    ,z+ h   *kz3);
ydot= y+ (h/6)*(ky1+ 2*ky2 + 2*ky3 + ky4);
zdot= z+ (h/6)*(kz1+ 2*kz2 + 2*kz3 + kz4);
%--------------------------------------------------------------------------
MR_states=[zdot ydot];
%--------------------------------------------------------------------------
end

