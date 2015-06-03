%parameter
m=0.05; %mass of the ball
e=1; %restit coeff
z0=1; %initial position
v0=-1; %initial velocity
g=9.81; %gravity
dt=0.000001; %time step
N=5000000; %number of step
R=160e-6; %radius of the ball
E=2.05e11; %Young modulus
nu=0.3; %Poisson's
Ep=1/(2*(1-nu^2)/E);
v=[v0];
z=[z0];
for nn=1:N
    vh=v(nn)-0.5*g*dt;
    z(nn+1)=z(nn)+vh*dt;
    if z(nn+1)>0
        a(nn+1)=-g;
        v(nn+1)=vh+0.5*a(nn+1)*dt;
    else
        a(nn+1)=-g+4/3*sqrt(R)*Ep*(-z(nn+1))^(3/2)/m;
        v(nn+1)=vh+0.5*a(nn+1)*dt;
    end
end
t=0:dt:N*dt;
plot(t,z,t,v);