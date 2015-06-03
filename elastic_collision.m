%parameter
m=0.05; %mass of the ball
e=1; %restit coeff
z0=1; %initial position
v0=-1; %initial velocity
g=9.81; %gravity
dt=0.05; %time step
N=1000; %number of step
v=[v0];
z=[z0];
for nn=1:N
    vh=v(nn)-0.5*g*dt;
    z(nn+1)=z(nn)+vh*dt;
    if z(nn+1)>0
        v(nn+1)=vh-0.5*g*dt;
    else
        v(nn+1)=-e*(vh-0.5*g*dt);
    end
end

%% Plot
t=0:dt:N*dt;
[hAx,~,~]=plotyy(t,z,t,v);

xlabel('Time (s)')
ylabel(hAx(1),'Position (m)') % left y-axis
ylabel(hAx(2),'Velocity (m/s)') % right y-axis