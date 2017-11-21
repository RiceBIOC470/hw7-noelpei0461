function gw2
dt=0.01;
interval=[0 10];
nstep=(interval(2)-interval(1))/dt;
xf2=0;
A=[];
for i=0.01:0.01:4 %take an a value from 0.01 to 4
    Xf=[];
    for ii=1:200 %random pick up 200 points
        x0=rand(1,1);
        sollx(1)=x0; %random initial value x0
                for ii=2:nstep
                    sollx(ii)=(sollx(ii-1)*(1-sollx(ii-1)))*a;
                    if xf2<sollx(ii)
                        xf2=sollx(ii);
                    end
                end
           Xf=[Xf,xf2];
    end
    A=[A;Xf];%append final Xf to A
end
a=0.01:0.01:4;
for iii=1:200
    plot(a,A(:,iii),'b.','MarkerSize',3) %plot Xf vs. a
    hold on;
end
end
