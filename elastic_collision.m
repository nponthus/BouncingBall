%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calcul Bouncing Ball elastic collision
%Nicolas
%03 06 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%parameter
m=0.05; %mass of the ball
e=1; %restit coeff
z0=1; %initial position
v0=-1; %initial velocity
g=9.81; %gravity
dt=0.05; %time step
N=10000; %number of step

%Velocity Verlet solving
v=[v0];
z=[z0];
Ec(1)=1/2*m*(v(1))^2;
Ep(1)=m*g*z(1);
Em(1)=Ec(1)+Ep(1);

for nn=1:N
    vh=v(nn)-0.5*g*dt;
    z(nn+1)=z(nn)+vh*dt;
    if z(nn+1)>0
        v(nn+1)=vh-0.5*g*dt;
    else
        v(nn+1)=-e*(vh-0.5*g*dt);
    end
    Ec(nn+1)=1/2*m*(v(nn+1))^2;
    Ep(nn+1)=m*g*z(nn+1);
    Em(nn+1)=Ec(nn+1)+Ep(nn+1);
end

%% Plot
t=0:dt:N*dt;
[hAx,~,~]=plotyy(t,z,t,v);

xlabel('Time (s)')
ylabel(hAx(1),'Position (m)') % left y-axis
ylabel(hAx(2),'Velocity (m/s)') % right y-axis

figure(2)
plot(t,Ep,t,Ec,t,Em)
% plot(t,Em)
title('Energy')
xlabel('Time (s)')
ylabel('Energy (J)')