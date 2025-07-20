%% Integracion de viga empotrada
close all; clear all; clc;

%% ===== Parametros fisicos viga unitaria =====
L = 1;                  % Longitud de la viga [m]
E = 0.5;
I = 0.5;
A = 1;                  % Area
EI = E*I;                 % Rigidez
M = 1;                  % Mass [kg]
rho = 200;
mu = 0;               % Amortiguamiento
alpha = 0.001;
beta = 0.0001;
P = 0.05;



K = (3*EI)/(L^3);

%% Parametros metodo numerico

Tiempo = 10000;
delta_t = 0.01;

N_step = round(Tiempo/delta_t);




%% nodes

ngl_node = 3;
n_nodes = 10;
n_elements = n_nodes - 1;

L_element = L/n_elements;

pos_nodes = zeros(n_nodes, 3);

N_gdl = n_nodes*ngl_node;

for i = 1:n_nodes

    pos_nodes(i, :) = [0, (i - 1)*L_element, 0];

end


%% Calculus of rigidity and mass matrix
K_global = zeros(n_nodes*ngl_node);
M_global = zeros(n_nodes*ngl_node);
gdle = zeros(n_elements, 6);

for i = 1:n_elements
    
    K_e = zeros(6);
    M_e = zeros(6);
    T_e = zeros(6);
    
    
    [K_e, Te, L_e] = rigidez_ele_flexion (E, A, I, pos_nodes(i, :), pos_nodes(i+1, :));
    [M_e] = mass_beam (A, rho, pos_nodes(i, :), pos_nodes(i+1, :));

    gdle(i, :) = [1, 2, 3, 4, 5, 6] + (i - 1)*[3, 3, 3, 3, 3, 3];

    idx = gdle(i, :);
    K_global(idx, idx) = K_global(idx, idx) + K_e;
    M_global(idx, idx) = M_global(idx, idx) + M_e;

end

C_global = alpha * M_global + beta * K_global;

%% Mass matrix

%% Beam empotrated

gdlR = [1, 2, 3];

gdlL = zeros(n_nodes*ngl_node - 3, 1);

for i = 1:(n_nodes*ngl_node - 3)
    gdlL(i) = i + 3;
end


%% 1. Condensacion estatica de matrices

KLL = K_global(gdlL, gdlL);            % Matriz de rigidez de gdl libres o reducida
KRR = K_global(gdlR, gdlR);            % Matriz de rigidez en gdl restringidos
KLR = K_global(gdlL, gdlR);            % Matriz de acoplamiento R-L
KRL = KLR';                     % Matriz de acoplamiento L-R


u_R = [0; 0; 0];              % Desplazamientos impuestos SPC
f_L = zeros(length(gdlL), 1);

f_L([(N_gdl - 5), (N_gdl - 4)]) = P*[1, 1];      % Vector de fuerzas externas en grados libres
u_L = KLL\f_L;

u_L_gs = Gauss_Seidel (KLL, f_L);

f_R = KRL*u_L + KRR*u_R;

U = zeros(N_gdl, 1);                       % Vector completo de desplazamientos

U(gdlL) = u_L;
U(gdlR) = u_R;

F = zeros(10, 1);                       % Vector completo de fuerzas

F(gdlL) = f_L;
F(gdlR) = f_R;

err_fint_fext = K*U -  F;

pos_plot = pos_nodes;


scale = 0.4/max(abs(U))

for i = 1:n_nodes

    ux = U(1 + 3*(i-1));
    uy = U(2 + 3*(i-1));
    utheeta = U(3 + 3*(i-1));
    pos_plot(i, :) = pos_plot(i, :) + scale*[ux, uy, utheeta];
end


%% Representacion



figure(1)
axis equal
hold on

%nodes init
plot(pos_nodes(:, 1), pos_nodes(:, 2), 'k')
plot(pos_nodes(:, 1), pos_nodes(:, 2), 'ko')
        

%nodes next
plot(pos_plot(:, 1), pos_plot(:, 2), 'r')
plot(pos_plot(:, 1), pos_plot(:, 2), 'ro')


%% simulation

U_t = U;
vel = zeros(length(gdlL), 1);
acc = zeros(length(gdlL), 1);

M_LL = M_global(gdlL, gdlL);
K_LL = K_global(gdlL, gdlL);
C_LL = C_global(gdlL, gdlL);
F_L = f_L;

% Figure 2: animación
figure(2)
axis([-1.3 1.3 -0.1 1.3*L])
hold on
title('Animation of the Deformed Beam. Scale Factor: 4')
xlabel('x')
ylabel('y')

% Posición inicial deformada
pos_dyn = pos_nodes;
for i = 1:n_nodes
    ux = U_t(1 + 3*(i-1));
    uy = U_t(2 + 3*(i-1));
    utheeta = U_t(3 + 3*(i-1));
    pos_dyn(i, :) = pos_dyn(i, :) + scale*[ux, uy, utheeta];
end

% Graficar una sola vez
curve = plot(pos_dyn(:,1), pos_dyn(:,2), 'b.-', ...
    'LineWidth', 2);
nodos_plot = plot(pos_dyn(:,1), pos_dyn(:,2), 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
drawnow


for step = 1:N_step

    acc = - M_LL \ (C_LL*vel + K_LL * u_L);
    vel = vel + acc * delta_t;
    u_L = u_L + vel * delta_t;

    % Reconstruir U completo
    U_t = zeros(N_gdl, 1);
    U_t(gdlL) = u_L;
    U_t(gdlR) = u_R;

    % Actualizar posiciones
    for i = 1:n_nodes
        ux = U_t(1 + 3*(i-1));
        uy = U_t(2 + 3*(i-1));
        utheeta = U_t(3 + 3*(i-1));
        pos_dyn(i, :) = pos_nodes(i, :) + scale*[ux, uy, utheeta];
    end

    % Actualizar curva de la animación
    curve.XData = pos_dyn(:,1);
    curve.YData = pos_dyn(:,2);

    nodos_plot.XData = pos_dyn(:,1);
    nodos_plot.YData = pos_dyn(:,2);
    drawnow limitrate
end
