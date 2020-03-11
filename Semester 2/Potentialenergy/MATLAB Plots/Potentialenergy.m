clear all
close all
figure ;
phic=csvread('potentialenergycorrection.csv');
phist=csvread('potentialenergyselftot.csv');
phisc=csvread('potentialenergyselfcorrect.csv');

phitc=csvread('potentialenergytotalcorrect.csv');
phitrap=csvread('potentialenergytrap.csv');



plot(phic,'.')
hold on 

plot(phist,'.')
hold on 

plot(phisc,'.')
hold on 
plot(phitc,'.')
hold on 
plot(phitrap,'.')

legend({'y = U self correction ','y = U self total','y = U self correct', 'y = UW +US correct','y = U Well total'})
title ('Potential energy per particles')
xlabel('Number of Particles') 
ylabel('Potential energy') 

hold on 

