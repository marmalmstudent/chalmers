# gamma0 = linspace(-0.19, 0.01, 100);
# gamma1 = linspace(-0.19, 0.01, 100);
# [g0, g1] = meshgrid(gamma0, gamma1);
# theta_m = 0.78063845074624338;
# z0 = ((1+g0)./(1-g0)).^2.*((1+g1)./(1-g1)).^2 - 1/2;
# # z1 = g0+g1+0.24748737341529153;
# z1 = g0+g1+0.16666666666666674;
# z2 = 2*g0 + 0.05/cos(theta_m)^3;
# z3 = 2*g1 + 3*0.05*(1/cos(theta_m)^3-1/cos(theta_m));
# hold on
# surf(g0,g1,abs(z0-z1))
# plot3(-0.069712888020348182, -0.10357390711963811, 0, 'o')
# surf(g0,g1,-abs(z2))
# surf(g0,g1,-abs(z3))

gamma0 = linspace(-0.15, -0.05, 100);
gamma1 = linspace(-0.25, -0.15, 100);
[g0, g1] = meshgrid(gamma0, gamma1);
theta_m = 0.7853981633974483;#0.809788457925;
gamma_m = 0.07731433;#0.07;
z0 = ((1+g0)./(1-g0)).^2.*((1+g1)./(1-g1)).^2 - 1/3.4;
z1 = 2*g0+2*g1+gamma_m*chebyshevpoly(1, 3, 1/cos(theta_m));
z2 = 2*g0 + gamma_m/cos(theta_m)^3;
z3 = 2*g1 + 3*gamma_m*(1/cos(theta_m)^3-1/cos(theta_m));
z4 = g0-g1/(3*sin(theta_m));
# data1 = z0.*(abs(z0) <= 1) + double(z0 > 1) - double(z0 < -1);
# data2 = z1.*(abs(z1) <= 1) + double(z1 > 1) - double(z1 < -1);
figure()
hold on
surf(g0,g1,z0)
# surf(g0,g1,abs(z4))
surf(g0,g1,z1)
v0 = -gamma_m/(2*cos(theta_m)^3);
v1 = -3/2*gamma_m*(1/cos(theta_m)^3-1/cos(theta_m));
f1 = ((1+v0)./(1-v0)).^2.*((1+v1)./(1-v1)).^2 - 1/3.4;
plot3(v0, v1, f1, 'ro')
plot3(-0.11785301, -0.18539951, 'bo')
surf(g0,g1,-abs(z2)-abs(z3))
daspect([1, 1, 1])
axis([-0.15, -0.05, -0.25, -0.15, -0.05, 0.05])
# surf(g0,g1,abs(z2+z3))