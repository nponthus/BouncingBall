%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calcul Bouncing Ball Hertzian contact
%Nicolas
%03 06 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%parameter
m=0.05; %mass of the ball
e=1; %restit coeff
z0=0; %initial position
v0=0; %initial velocity
g=9.81; %gravity
dt=0.000001; %time step
df=0.05; %length of the simulation (m)
R=160e-6; %radius of the ball
E=2.05e11; %Young modulus
nu=0.3; %Poisson's ratio
Er=1/(2*(1-nu^2)/E); %reduced modulus
K=4/3*sqrt(R)*Er/m;
VV=[0.001,0.002,0.005,0.007,0.01,0.012,0.015,0.017,0.02,0.025,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1];
lmc=0.000150;
Rq=0.000010;
load('HC.mat');
for jj=1:length(VV)
    V=VV(jj);
    tf=df/V;
    t=0:dt:tf;
    N=length(t);
    X=V*(0:dt:(N+1)*dt);
    h1=interp1(x,h1f,X);
    eta=0.001;
    c=4*sqrt(K*m^2)*eta;
    p=3;
    ds=(g/K)^(2/3);
    
    %coeff passage base normalisee/base relative
    Cz=Rq;
    Cv=Rq*V*2*pi/lmc;
    Ct=lmc/(2*pi*V);

    %Velocity Verlet solving
    v=[v0];
    z=[z0];
    a=[-g];
    Ec(1)=1/2*m*(v(1))^2;
    Ep(1)=m*g*z(1);
    Ek(1)=0;
    Em(1)=Ek(1)+Ec(1)+Ep(1);

    for nn=1:N-1
        vh=v(nn)+0.5*a(nn)*dt;
        z(nn+1)=z(nn)+vh*dt;
        if (z(nn+1)-h1(nn+1))>0
            a(nn+1)=-g;
            Ek(nn+1)=0;
        else
            a(nn+1)=-g+K*(h1(nn+1)-z(nn+1))^(3/2)+c/m*((h1(nn+2)-h1(nn))/(2*dt)-vh)*((h1(nn+1)-z(nn+1))/ds)^(p/2);
            Ek(nn+1)=2/5*K*m*(h1(nn+1)-z(nn+1))^(5/2);
        end
        v(nn+1)=vh+0.5*a(nn+1)*dt;
        Ec(nn+1)=1/2*m*(v(nn+1))^2;
        Ep(nn+1)=m*g*z(nn+1);
        Em(nn+1)=Ek(nn+1)+Ec(nn+1)+Ep(nn+1);
    end
    tn=1/Ct*t;
    zn=1/Cz*z;
    vn=1/Cv*v;
save(strcat('mouvement',num2str(V*1000),'mms_norm.mat'),'tn','zn','vn')
end
% %% Plot
% t=0:dt:N*dt;
% figure(1)
% [hAx,~,~]=plotyy(t,z,t,v);
% xlabel('Time (s)')
% ylabel(hAx(1),'Position (m)') % left y-axis
% ylabel(hAx(2),'Velocity (m/s)') % right y-axis
% 
% figure(2)
% plot(t,Em,t,Ep,t,Ec,t,Ek)
% %plot(t,Em)
% title('Energy')
% xlabel('Time (s)')
% ylabel('Energy (J)')
% legend('Em','Ep','Ec','Ek')
% 
% figure(3)
% plot(X/V,h1,t,z)
% title('displacements')
% xlabel('time (s)')
% ylabel('position(m)')
% legend('surface','ball')