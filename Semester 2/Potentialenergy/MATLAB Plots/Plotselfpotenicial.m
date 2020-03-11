clear all
close all
figure ;
phi=csvread('Self.csv');
phi_new = reshape(phi, 201, 100);
h = surf(phi_new)
set(h,'LineStyle','none') 
title ('self Potential')
xlabel('Nr') 
ylabel('Nz') 
zlabel('Potential (V)') 

figure ;

phi11=csvread('correctionpotential11.csv');
hold on 
phi12=csvread('correctionpotential12.csv');

Nr=1:1:99;

plot(Nr,phi11)
hold on 
plot (Nr,phi12)
title ('self Potential')
xlabel('Nr') 
ylabel('Potentials') 
legend({'y = phi11','y = phi12'})
