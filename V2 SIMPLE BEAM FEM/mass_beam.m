function [M_e] = mass_beam(A, rho, n_1, n_2)
    % n_1, n_2: coordenadas de nodos del elemento

    deltax_e = n_2(1) - n_1(1);
    deltay_e = n_2(2) - n_1(2);
    L_e = sqrt(deltax_e^2 + deltay_e^2);

    m = rho * A * L_e;  % masa total del elemento

    % Matriz de masa consistente en coordenadas locales
    M_local = (m / 420) * [
        140  0     0     70   0     0;
        0    156   22*L_e 0    54   -13*L_e;
        0    22*L_e 4*L_e^2 0  13*L_e -3*L_e^2;
        70   0     0    140   0     0;
        0    54    13*L_e 0   156   -22*L_e;
        0   -13*L_e -3*L_e^2 0 -22*L_e 4*L_e^2
    ];

    % Matriz de transformación
    alfa = atan2(deltay_e, deltax_e);
    R = [cos(alfa), sin(alfa), 0; 
        -sin(alfa), cos(alfa), 0;
         0,         0,         1];
    T_e = blkdiag(R, R);  % matriz 6x6

    % Transformar a coordenadas globales
    M_e = T_e' * M_local * T_e;
end