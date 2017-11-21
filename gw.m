function dx=gw(x0,a)
dt=0.01;
interval=[0 10/a];
nstep=(interval(2)-interval(1))/dt;
sollx(1)=x0;
for ii=2:nstep
    sollx(ii)=sollx(ii-1)+(sollx(ii-1)*(1-sollx(ii-1)))*a*dt;
end
tt=linspace(interval(1),interval(2), nstep);
plot(tt,sollx)
t=[];
b=1;
for i=1:nstep
    if sollx(i)>=0.99*a/4
        t(b)=tt(i);
        b=b+1;
    end
end
min(t)
end
