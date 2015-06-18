function [ A , Ek] = acceleration(z_ballsc, rho_balls, theta_balls,v_balls,eta)
%Returns acceleration and elastic potential Energy
global m g K I J ds p
z_balls=reshape(z_ballsc,1,length(z_ballsc));
v_balls=reshape(v_balls,1,length(v_balls));
    ind=find(z_balls <=0);
    if isempty(ind)==false
        A(1)=-g+sum(K(ind).*(-z_balls(ind)).^(3/2))/m-sum(4*eta*sqrt(K(ind)*m/length(ind))/m.*v_balls(ind).*(-z_balls(ind)/ds(ind)).^(p/2));
        A(2)=sum(K(ind).*(-z_balls(ind)).^(3/2).*rho_balls(ind).*sin(theta_balls(ind)))/I-sum(4*eta*sqrt(K(ind)*(m/length(ind)))/I.*v_balls(ind).*(-z_balls(ind)/ds(ind)).^(p/2).*rho_balls(ind).*sin(theta_balls(ind)));
        A(3)=-sum(K(ind).*(-z_balls(ind)).^(3/2).*rho_balls(ind).*cos(theta_balls(ind)))/J+sum(4*eta*sqrt(K(ind)*(m/length(ind)))/J.*v_balls(ind).*(-z_balls(ind)/ds(ind)).^(p/2).*rho_balls(ind).*cos(theta_balls(ind)));
        Ek=sum(2/5*K(ind).*(-z_balls(ind)).^(5/2));
    else
        A(1)=-g;
        A(2)=0;
        A(3)=0;
        Ek=0;
    end
end