% ----- Parameters -----

% Number of grid points on the inflow boundary (the grid points in the 
% rest of the domain is scaled accordingly)
m = 15;

% Reynolds number
%Re = 80;
Re = 20:30:200;
% End time
Tend = 6;

% Length of outflow pipe
Lout = 2;

% Show animation and streamlines of the solution? True or false. 
animate = false;

% ----- Main function -----
x_ra = zeros(length(Re), 1);
for i = 1:length(Re)
    x_ra(i) = incomp(m,Re(i),Tend,Lout,animate);
end
%%

plot(Re, x_ra, '-o')

xlabel('Reynolds number')
ylabel('Reattachment length')
title('Reattachment length vs Reynolds number')
grid on
% ----- End of main function -----
%x_ra är kortaste distancen där hastigheten i x-led är helt positiv
%Re = linspace(20, 200, 30);
%För varje Re specificera spatial resolution, end time och outflow pipens
%längd.