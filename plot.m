clear;
close all;

load("Computation\data")

radius = 0.5;
Re = 200.0;
uInf = 1.0;

% Maillage espace
sizeR = 40;
sizeTheta = 40;

% Maillage en coordonnée polaire
lengthR = 19.0 * radius;
stepSizeR = lengthR / (sizeR - 1.0);
r = linspace(radius, radius + lengthR, sizeR);

lengthTheta = 2.0 * pi;
stepSizeTheta = lengthTheta / (sizeTheta - 1.0);
theta = linspace(0, lengthTheta, sizeTheta);

x = r' .* cos(theta);
y = r' .* sin(theta);

% APPROXIMATION DES VITESSES PAR DIFFERENCES FINIES
% CALCUL DE UR
ur = zeros(sizeR,sizeTheta);
for i=1:sizeR
    for j = 2:sizeTheta-1
        ur(i,j)  = ( data(i,j+1)-data(i,j-1) ) / (2*r(i)*stepSizeTheta); 
    end
    ur(i,1) = ( data(i,2) - data(i,1) )/ (r(i)*stepSizeTheta); 
    ur(i,sizeTheta) = ( data(i,sizeTheta) - data(i,sizeTheta-1))/ (r(i)*stepSizeTheta); 
end

% CALCUL DE UT
ut = zeros(sizeR,sizeTheta);
for j = 1:sizeR 
    for i=2:sizeTheta-1
       ut(i,j) = - ( data(i+1,j) - data(i-1,j) ) / (2*stepSizeR);  
    end
    ut(1,j) = - ( data(2,j) - data(1,j)   ) / stepSizeR;  
    ut(sizeR,j) = - ( data(sizeR,j) - data(sizeR-1,j) ) / stepSizeR;  
end

ur = ur';
ut = ut';

% PROJECTION DANS LA BASE CARTESIENNE
ux = ur.*cos(theta') - ut.*sin(theta');   
uy = ur.*sin(theta') + ut.*cos(theta');

%CALCUL DE LA NORME DE LA VITESSE
u = sqrt(ux .^ 2 + uy .^ 2);

% GRAPHIQUE DES VITESSES NUMERIQUES
fig = pcolor(x,y,u');
axis("equal")
set(fig, 'EdgeColor', 'none');