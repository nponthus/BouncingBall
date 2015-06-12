clear all
d=0.1;
V=0.001;
dt=0.000001;
tf=d/V;
t=0:dt:tf;
Ra=0.00001;
lamc=0.000150;
fc=V/lamc
h1(1)=0;
h1f(1)=0;
hp(1)=0;
a(1)=0;
z=0.25;
w=fc*2*pi;
for nn=1:(length(t)-1)
    hph=hp(nn)+1/2*a(nn)*dt;
    h1(nn+1)=h1(nn)+hph*dt;
    a(nn+1)=-2*z*w*hph-w^2*h1(nn+1)+randn*4*sqrt(z*w*w*w/dt);
    hp(nn+1)=hph+1/2*a(nn+1)*dt;
    h1f(nn+1)=h1f(nn)-1.3*w*h1f(nn)*dt+h1(nn+1)*dt;
end
h1=(h1-mean(h1))*Ra/std(h1);
h1f=(h1f-mean(h1f))*Ra/std(h1f);
plot(V*t,h1,V*t,h1f)
%dfittool(h1f)
x=V*t;
save('H1.mat','h1f','x')
%run('J:\données TFE\profilo\rayon.m')