%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calcul Bouncing Ball 3 balls on a plate
%Nicolas
%03 06 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameter
m=0.05; %mass of the ball
I=m/12*(0.025^2+0.01^2);
J=m/12*(0.025^2+0.01^2);
zg0=0.02; %initial position of the center of mass
vg0=0; %initial velocity of the center of mass
phi0=0.3; %initial angle along x
phip0=0; %initial rotational velocity along x
psi0=0.4; %initial angle along y
psip0=0; %initial rotational velocity along y
g=9.81; %gravity
dt=0.000001; %time step
N=1000000; %number of step
R=[160e-6 160e-6 160e-6]; %radius of the balls
E=2.05e11; %Young modulus
nu=0.3; %Poisson's ratio
Ep=1/(2*(1-nu^2)/E); %reduced modulus
K=4/3*sqrt(R)*Ep; %Hertz force constant
rho=[0.01 0.01 0.01]; %localisation of the contact point, for now an equilateral triangle
theta=[0 2*pi/3 4*pi/3];

%Initialisation
zg=[zg0]; %G position
vg=[vg0]; %G velocity
phi=[phi0]; %X rotation
phip=[phip0]; %X rotational velocity
psi=[psi0]; %Y rotation
psip=[psip0]; %Y rotational velocity
z=zeros(3,N+1); % table of the position of contact point
z(:,1)=zg0+phi0*rho.*sin(theta)-psi0*rho.*cos(theta); %initial position of the contact point supposed to be positiiv

%Velocity Verlet solving
for nn=1:N
    vh(1)=vg(nn)-1/2*g*dt;
    vh(2)=phip(nn);
    vh(3)=psip(nn);
    zg(nn+1)=zg(nn)+vh(1)*dt;
    phi(nn+1)=phi(nn)+vh(2)*dt;
    psi(nn+1)=psi(nn)+vh(3)*dt;
    z(:,nn+1)=zg(nn+1)+phi(nn+1)*rho.*sin(theta)-psi(nn+1)*rho.*cos(theta);
    ind=find(z(:,nn+1)<=0);
    if isempty(ind)==false
        a(1)=-g+sum(K(ind)*(-z(ind,nn+1)).^(3/2))/m;
        a(2)=sum(K(ind)*(-z(ind,nn+1)).^(3/2).*rho(ind).*sin(theta(ind)))/I;
        a(3)=-sum(K(ind)*(-z(ind,nn+1)).^(3/2).*rho(ind).*cos(theta(ind)))/J;
    else
        a(1)=-g;
        a(2)=0;
        a(3)=0;
    end
    vg(nn+1)=vh(1)+1/2*a(1)*dt;
    phip(nn+1)=vh(2)+1/2*a(2)*dt;
    psip(nn+1)=vh(3)+1/2*a(3)*dt;
end

%% Plot
t=0:dt:N*dt;
figure(1)
[hAx,~,~]=plotyy(t,zg,t,vg);
title('movement of the center of mass')
xlabel('Time (s)')
ylabel(hAx(1),'Position (m)') % left y-axis
ylabel(hAx(2),'Velocity (m/s)') % right y-axis

figure(2)
plot(t,z(1,:),t,z(2,:),t,z(3,:));
title('position of the contact point')
xlabel('Time (s)')
ylabel('Position (m)')
