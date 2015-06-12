function [ F,Fp ] = static( x,A,t,r,h )
%static equation of three balls
global m g
x=reshape(x,length(x),1);
A=reshape(A,length(A),1);
t=reshape(t,length(t),1);
r=reshape(r,length(r),1);
h=reshape(h,length(h),1);
F(1,1)=sum(A.*(h-x).^(3/2))-m*g;
F(2,1)=sum(A.*r.*sin(t).*(h-x).^(3/2));
F(3,1)=sum(A.*r.*cos(t).*(h-x).^(3/2));
Fp(1,:)=A.*(-3/2*(h-x).^(1/2));
Fp(2,:)=A.*r.*sin(t).*(-3/2*(h-x).^(1/2));
Fp(3,:)=A.*r.*cos(t).*(-3/2*(h-x).^(1/2));
end
