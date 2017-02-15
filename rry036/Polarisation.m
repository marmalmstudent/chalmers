T=1;
t=linspace(0,T*3/4,100);
w=2*pi/T;
X=sqrt(3)*exp(1i*w*t);
Y=-exp(1i*w*t);
plot(real(X), real(Y))
axis([-2, 2, -2, 2])