%HW7
%GB comments
1a 100
1b 100
1c 100
1d 50 The wrong plot is generated and the comments do not explain the behavior accurately.  Your function should have returned a bifurcation diagram that illustrates the emergence of new fixed points as a function of ‘a.’
2a. 50 wrong equations used. It should be [V/(1+x(2)^4)-x(1); V/(1+x(1)^4)-x(2)];
2b. 100 correctly generated plots, but the equations used are wrong. I will give full credit
2c  100 same as 2b. 
overall: 86


% Problem 1: Modeling population growth
% The simplest model for a growing population assumes that each current
% individual has equal likelihood to divide, which yields a differential
% equation dx/dt = a*x where a is the division rate. It is easy to see that
% this will grow exponentially without bound. A simple modification is to
% assume that the growth rate slows done as the population reaches some
% maximum value N so that dx/dt = a*x*(1-x/N). Defining X = x/N, we have 
% dX/dt = a*X*(1-X).  
% Part 1. This equation has two fixed points at 0 and 1. Explain the
% meaning of these two points.

% The fix point at 0 means when the initial population is 0, the total
% population will not change according to time. Because in this case, no
% initial population is shown and thus it can not grow. So it is a fixed
% point.

% The time point at 1 means the population growth reaches its maximum capacity N. The
% population is in a saturated stage and cannot grow any more. Thus the it
% is another fixed point.

% Part 2: Evaluate the stability of these fixed points. Does it depend on
% the value of the parameter a? 

a=2;
dt=0.01;
interval=[0 10/a];
nstep=(interval(2)-interval(1))/dt;
sollx(1)=0;
for ii=2:nstep
    sollx(ii)=sollx(ii-1)+(sollx(ii-1)*(1-sollx(ii-1)))*a*dt;
end
tt=linspace(interval(1),interval(2), nstep);
plot(tt,sollx)

a=5;
dt=0.01;
interval=[0 10/a];
nstep=(interval(2)-interval(1))/dt;
sollx(1)=0;
for ii=2:nstep
    sollx(ii)=sollx(ii-1)+(sollx(ii-1)*(1-sollx(ii-1)))*a*dt;
end
tt=linspace(interval(1),interval(2), nstep);
plot(tt,sollx)

a=10;
dt=0.01;
interval=[0 10/a];
nstep=(interval(2)-interval(1))/dt;
sollx(1)=0;
for ii=2:nstep
    sollx(ii)=sollx(ii-1)+(sollx(ii-1)*(1-sollx(ii-1)))*a*dt;
end
tt=linspace(interval(1),interval(2), nstep);
plot(tt,sollx)

% No, the stability of these fixed points are very good and independent
% from a. These are the initial and final stage of the growth. So,
% regardless of the given growth factor a, these fixed point will not
% change. However, the slope of does depend on the a. The higher the a
% value, the faster the population reaches final stage.

% Part 3: Write a function that takes two inputs - the initial condition x0
% and the a parameter and integrates the equation forward in time. Make
% your code return two variables - the timecourse of X and the time
% required for the population to reach 99% of its maximum value.
x0=0.5;
a=2;
gw(x0,a) % gw is the function I generated for this problem.


% Part 4: Another possible model is to consider discrete generations
% instead allowing the population to vary continuously. e.g. X(t+1) = a*
% X(t)*(1-X(t)). Consider this model and vary the a parameter in the range 0
% < a <= 4. For each value of a choose 200 random starting points  0 < x0 < 1 
% and iterate the equation forward to steady state. For each final
% value Xf, plot the point in the plane (a,Xf) so that at the end you will
% have produced a bifucation diagram showing all possible final values of
% Xf at each value of a. Explain your results. 

gw2; 
% The plot shows that as the a increases, the final value Xf also
% increases and is proportional to a. The result corresponds to the
% expectation, where the final value is a/4, since the maximum value of
% (1-x)*x is 0.25.


% Problem 2. Genetic toggle switches. 
% Consider a genetic system of two genes A and B in which each gene
% product represses the expression of the other. Make the following
% assumptions: 
% a. Repression is cooperative:  each promotor region of one gene has 4
% binding sites for the other protein and that all of these need to be
% occupied for the gene to be repressed. 
% b. You can ignore the intermediate mRNA production so that the product of
% the synthesis of one gene can be assumed to directly repress the other
% c. the system is prefectly symmetric so that the degradation
% times, binding strengths etc are the same for both genes. 
% d. You can choose time and concentration scales so that all Michaelis
% binding constants and degradation times are equal to 1. 
%
% Part 1. Write down a two equation model (one for each gene product) for
% this system. Your model should have one free parameter corresponding to the
% maximum rate of expression of the gene, call it V. 
% dX1/dt=(V+(V-k)X2^4)/(1 + X2^4)-X2;dX2/dt=(V+(V-k)X1^4)/(1 + X1^4)-X1;



% Part 2. Write code to integrate your model in time and plot the results for V = 5 for two cases, 
% one in which A0 > B0 and one in which B0 > A0. 
% A0>B0
k=1; 
V=5; 
A0=10; 
B0=2;
dxA=@(t,B) (V+(V-k).*B^4)/(1 + B^4)-B;
dxB=@(t,A) (V+(V-k).*A^4)/(1 + A^4)-A;
solA1=ode23(dxA, [0 10], B0);
solB1=ode23(dxB, [0 10], A0);
figure; hold on;
plot(solA1.x,solA1.y,'r-'); 
hold on;
plot(solB1.x,solB1.y,'g-');
hold off;
xlabel('time'); ylabel('expression');
legend('A', 'B');
% A0<B0
k=1; 
V=5; 
A0=2; 
B0=10;
dxA=@(t,B) (V+(V-k).*B^4)/(1 + B^4)-B;
dxB=@(t,A) (V+(V-k).*A^4)/(1 + A^4)-A;
solA1=ode23(dxA, [0 10], B0);
solB1=ode23(dxB, [0 10], A0);
figure; hold on;
plot(solA1.x,solA1.y,'r-'); 
hold on;
plot(solB1.x,solB1.y,'g-'); 
hold off;
xlabel('time'); ylabel('expression');
legend('A', 'B');
% Part 3. By any means you want, write code to produce a bifurcation diagram showing all
% fixed points of the system as a function of the V parameter.

% dA/dt = 0 = -B^5+(V-k)*B^4+0*B^3+0*B^2-B+V
figure(1); 
hold on;
k=2; 
for V=0:0.05:5
    polycoeff=[-1 (V-k) 0 0 -1 V]
    rts=roots(polycoeff);
    rts=rts(imag(rts)==0);
    plot(V*ones(length(rts),1),rts,'r.');
end
hold off;
xlabel('V value'); ylabel('fixed points');
