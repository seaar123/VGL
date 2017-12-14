function [xdot,df_dx,df_du]=Dynamic_MRobot(t,x,u)
tau=u';
%--------------------------------------------------------------------------
%Robot States
xc=x(1);%x position
yc=x(2);%y position
th=x(3);%azmith angle
phR=x(4);%Right angular position
phL=x(5);%left angular position
v=x(6);% linear Velocity
w=x(7);% angular Velocity
V=[v;w];
%--------------------------------------------------------------------------
%Robot Parameters
mT=x(8);
mW=x(9);
r=x(10);
b=x(11);
d=x(12);
Iyy=x(13);
IT=x(14);
% ignore friction *0
Fv=x(15); % *0
Fd=x(16); % *0
%--------------------------------------------------------------------------
%Robot Kinematics
S=[cos(th) -d*sin(th);
   sin(th)  d*cos(th);
   0        1;
   1/r      b/r;
   1/r     -b/r];
%--------------------------------------------------------------------------
% next states ( position )
xdot(1:5)=S*V;
phRdot=xdot(4);
phLdot=xdot(5);
%--------------------------------------------------------------------------
%dynamics
%------------
m11=mT+2*Iyy/(r^2);
m22=mT*d^2+IT+2*Iyy*(b^2)/(r^2)-4*mW*d^2;
M=[m11  0 ;
   0  m22];
%------------
Vm=[0  -d*w*(mT-2*mW);
    d*w*(mT-2*mW)   0];
%------------ friction
F=(1/r)*[Fv*(phRdot+phLdot)+Fd*(sign(phRdot)+sign(phLdot));
      b*(Fv*(phRdot-phLdot)+Fd*(sign(phRdot)-sign(phLdot)))];
%------------gravity
G=[0;0];
%------------desired torque vescuing
td=[0;0];
%------------
B=[1/r  1/r;
   b/r -b/r];
%--------------------------------------------------------------------------
% next states ( velocity )
xdot(6:7)=M^(-1)*(B*tau-Vm*V-td-G-F);
%--------------------------------------------------------------------------
Im11=1/m11;
Im22=1/m22;
df_dx=[0 -2*w*Im11*d*(mT-2*mW);w*Im22*d*(mT-2*mW) v*Im22*d*(mT-2*mW)];
df_du=M^(-1)*B; % becuase it is affine system ( xdot=f(x)+g(x)u )

