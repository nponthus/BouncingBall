%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calcul Bouncing Ball 3 balls on a plane with damping
%Nicolas
%03 06 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
tic
global m g K I J ds p
%parameter
m=0.05; %mass of the ball
I=m/12*(0.025^2+0.01^2);
J=m/12*(0.025^2+0.01^2);
zg0=0.00001; %initial position of the center of mass
vg0=0; %initial velocity of the center of mass
phi0=0; %initial angle along x
phip0=0; %initial rotational velocity along x
psi0=0; %initial angle along y
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
Mbase=[1 rho(1)*sin(theta(1)) rho(1)*cos(theta(1));1 rho(2)*sin(theta(2)) rho(2)*cos(theta(2));1 rho(3)*sin(theta(3)) rho(3)*cos(theta(3))];
eta=0.001;
p=0;
ds=(m/3*g./K).^(2/3) 

%Initialisation
zg=[zg0]; %G position
vg=[vg0]; %G velocity
phi=[phi0]; %X rotation
phip=[phip0]; %X rotational velocity
psi=[psi0]; %Y rotation
psip=[psip0]; %Y rotational velocity
z=zeros(3,N+1); % table of the position of contact point
z(:,1)=zg0+phi0*rho.*sin(theta)-psi0*rho.*cos(theta); %initial position of the contact

% test=find(z(:,1)<0);
% if isempty(test)==false
%     error('wrong initial position')
% end

[a,Ek(1)]=acceleration(z(:,1),rho,theta);
Ec(1)=1/2*m*vg(1)^2+1/2*I*phip(1)^2+1/2*J*psip(1)^2;
Ep(1)=m*g*zg(1);
Em(1)=Ec(1)+Ek(1)+Ep(1);

%Velocity Verlet solving
for nn=1:N
    vh(1)=vg(nn)+1/2*a(1)*dt;
    vh(2)=phip(nn)+1/2*a(2)*dt;
    vh(3)=psip(nn)+1/2*a(3)*dt;
    zg(nn+1)=zg(nn)+vh(1)*dt;
    phi(nn+1)=phi(nn)+vh(2)*dt;
    psi(nn+1)=psi(nn)+vh(3)*dt;
    z(:,nn+1)=zg(nn+1)+phi(nn+1)*rho.*sin(theta)-psi(nn+1)*rho.*cos(theta);
    vhballs=Mbase*vh';
    
    [a,Ek(nn+1)]=dampedacc(z(:,nn+1),rho,theta,vhballs,eta);
    
    vg(nn+1)=vh(1)+1/2*a(1)*dt;
    phip(nn+1)=vh(2)+1/2*a(2)*dt;
    psip(nn+1)=vh(3)+1/2*a(3)*dt;
    Ec(nn+1)=1/2*m*(vg(nn+1))^2+1/2*I*(phip(nn+1))^2+1/2*J*(psip(nn+1))^2;
    Ep(nn+1)=m*g*zg(nn+1);
    Em(nn+1)=Ec(nn+1)+Ek(nn+1)+Ep(nn+1);
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

figure(3)
plot(t,phi,t,psi);
title('Rotation of the slider')
xlabel('Time (s)')
ylabel('Rotation (rad)')

figure(4)
%plot(t,Em,t,Ep,t,Ec,t,Ek)
plot(t,Em)
title('Energy')
xlabel('Time (s)')
ylabel('Energy (J)')
legend('Em','Ep','Ec','Ek')
toc