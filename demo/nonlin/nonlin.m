nu = 1000;
u = linspace(0,2*pi,nu)';
xtest = linspace(-2,2,400)';
nx = length(xtest);

% Passband function
%fp = @(x) 2./(1 + exp(-x))-1;
fp = @(x) x - 0.1*x.^3;

% Baseband equivalent function
cosu = cos(u)';
fb = 2*mean(cosu.*fp(xtest*cosu),2);

plot(xtest, [fp(xtest) fb]);


