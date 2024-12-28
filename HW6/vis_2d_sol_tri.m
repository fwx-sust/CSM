clear all; clc; clf;

load("HEAT.mat");

hh_x = 2.0 / n_el_x;
hh_y = 2.0 / n_el_y;

n_np_x = 0.5*n_el_x + 1;
n_np_y = 0.5*n_el_y + 1;

[X, Y] = meshgrid(0 : hh_x : 1, 0 : hh_y : 1);
Z = reshape(disp, n_np_x, n_np_y)';
surf(X, Y, Z);

shading interp

az = -61;
el = 20;
view(az,el);

% EOF