

phi=csvread('hered0.csv');
phi_new = reshape(phi, 401, 150);

figure;
    surf(phi_new)
    shading interp,colormap(jet(1024)),colorbar
    title('plot of the solution of the poisson boltzmann equation')
    zlabel('chargedensity')
    xlabel('r')
    ylabel('z')
    
   
    phi=csvread('hered1.csv');
phi_new = reshape(phi, 401, 150);

figure;
    surf(phi_new)
    shading interp,colormap(jet(1024)),colorbar
    title('plot of the solution of the poisson boltzmann equation')
    zlabel('chargedensity')
    xlabel('r')
    ylabel('z')
    
  % fileID = fopen('herebin.bin');
  % A = fread(fileID)
  phi=csvread('hered2.csv');
phi_new = reshape(phi, 401, 150);

figure;
    surf(phi_new)
    shading interp,colormap(jet(1024)),colorbar
    title('plot of the solution of the poisson boltzmann equation')
    zlabel('chargedensity')
    xlabel('r')
    ylabel('z')
    