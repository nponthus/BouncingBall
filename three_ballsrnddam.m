%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calcul Bouncing Ball 3 balls on a plane
%Nicolas
%03 06 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
tic
global m g K I J ds p
%parameter
m=0.05; %mass of the ball
I=m/12*(0.025^2+0.01^2);
J=m/12*(0.025^2+0.01^2);
vg0=0; %initial velocity of the center of mass
phip0=0; %initial rotational velocity along x
psip0=0; %initial rotational velocity along y
g=9.81; %gravity
dt=0.000001; %time step
xf=0.005;
R=[80e-6 80e-6 80e-6]; %radius of the balls
E=2.05e11; %Young modulus
nu=0.3; %Poisson's ratio
Er=1/(2*(1-nu^2)/E); %reduced modulus
K=4/3*sqrt(R)*Er; %Hertz force constant
rho=[0.01 0.01 0.01]; %localisation of the contact point, for now an equilateral triangle
theta=[0 2*pi/3 4*pi/3];
Mbase=[1 rho(1)*sin(theta(1)) rho(1)*cos(theta(1));1 rho(2)*sin(theta(2)) rho(2)*cos(theta(2));1 rho(3)*sin(theta(3)) rho(3)*cos(theta(3))];
eta=0.035;
ds=(m/3*g./K).^(2/3);
p=3;
load H1.mat;
load H2.mat;
load H3.mat;

V=0.007
tf=xf/V;
t=0:dt:tf;
X=V*t;
N=length(t)-1;

h(1,:)=interp1(x,h1f,X);
h(2,:)=interp1(x,h2f,X);
h(3,:)=interp1(x,h3f,X);
hp=(h(:,2:end)-h(:,1:end-1))./dt;

%Initialisation
vg=[vg0]; %G velocity
phip=[phip0]; %X rotational velocity
psip=[psip0]; %Y rotational velocity
z=zeros(3,N+1); % table of the position of contact point
z(:,1)=equilibre(@static,h(:,1)-0.1*[1; 1; 1],K,theta,rho,h(:,1),1e-5,100);
Mbase=[1 rho(1)*sin(theta(1)) rho(1)*cos(theta(1));1 rho(2)*sin(theta(2)) rho(2)*cos(theta(2));1 rho(3)*sin(theta(3)) rho(3)*cos(theta(3))];
p0=Mbase^(-1)*z(:,1);
zg=[p0(1)];
phi=[p0(2)];
psi=[p0(3)];

% test=find(z(:,1)<0);
% if isempty(test)==false
%     error('wrong initial position')
% end

[a,Ek(1)]=accelerationrnd(z(:,1),rho,theta,h(:,1));
Ec(1)=1/2*m*vg(1)^2+1/2*I*phip(1)^2+1/2*J*psip(1)^2;
Ep(1)=m*g*zg(1);
Em(1)=Ec(1)+Ek(1)+Ep(1);
ag=[a(1)];

%Velocity Verlet solving
for nn=1:N
    vh(1)=vg(nn)+1/2*a(1)*dt;
    vh(2)=phip(nn)+1/2*a(2)*dt;
    vh(3)=psip(nn)+1/2*a(3)*dt;
    zg(nn+1)=zg(nn)+vh(1)*dt;
    phi(nn+1)=phi(nn)+vh(2)*dt;
    psi(nn+1)=psi(nn)+vh(3)*dt;
    z(:,nn+1)=zg(nn+1)+phi(nn+1)*rho.*sin(theta)-psi(nn+1)*rho.*cos(theta);
    delta=z(:,nn+1)-h(:,nn+1);
    vr=Mbase*vh'-hp(:,nn);

    [a,Ek(nn+1)]=dampedacc(delta,rho,theta,vr,eta);

    ag(nn+1)=a(1);
    vg(nn+1)=vh(1)+1/2*a(1)*dt;
    phip(nn+1)=vh(2)+1/2*a(2)*dt;
    psip(nn+1)=vh(3)+1/2*a(3)*dt;
    Ec(nn+1)=1/2*m*(vg(nn+1))^2+1/2*I*(phip(nn+1))^2+1/2*J*(psip(nn+1))^2;
    Ep(nn+1)=m*g*zg(nn+1);
    Em(nn+1)=Ec(nn+1)+Ek(nn+1)+Ep(nn+1);
end


%% Plot
% figure(1)
% [hAx,~,~]=plotyy(t,zg,t,vg);
% title('movement of the center of mass')
% xlabel('Time (s)')
% ylabel(hAx(1),'Position (m)') % left y-axis
% ylabel(hAx(2),'Velocity (m/s)') % right y-axis
% 
% figure(2)
% plot(t,z(1,:),t,z(2,:),t,z(3,:));
% title('position of the contact point')
% xlabel('Time (s)')
% ylabel('Position (m)')
% 
% figure(3)
% plot(t,phi,t,psi);
% title('Rotation of the slider')
% xlabel('Time (s)')
% ylabel('Rotation (rad)')
% 
% figure(4)
% %plot(t,Em,t,Ep,t,Ec,t,Ek)
% plot(t,Em)
% title('Energy')
% xlabel('Time (s)')
% ylabel('Energy (J)')
% legend('Em','Ep','Ec','Ek')
% 
% figure(5)
% plot(t,zg*1e6,t,vg*1000,t,ag)
% title('movement of the center of mass')
% xlabel('Time (s)')
% legend('displacement(µm)','velocity (mm/s)','acceleration (m\^2/s)')
% toc