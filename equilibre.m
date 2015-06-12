function [ x,iter ] = equilibre( fun,x0,A,t,r,h,eps,imax )
%give the altitude of the balls at t=0
global m g
[f,fp]=fun(x0,A,t,r,h);
iter=1;
x=reshape(x0,length(x0),1);
while abs(f(1))>eps && iter<=imax
    x=x-fp^(-1)*f;
    [f,fp]=fun(x,A,t,r,h);
    iter=iter+1;
end
end

