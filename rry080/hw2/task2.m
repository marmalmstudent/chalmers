clear all;
load("constants.mat")  % the constants from the table
load("sarlab.mat")

L = (R0_near + 349*dR)*theta_a;
printf("Synthetic aperture: %f m = %f pixels\n",
	L, L/dA);