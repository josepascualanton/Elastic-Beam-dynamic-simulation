function [K_e, Te, L_e] = rigidez_ele_flexion (E, A, I, n_1, n_2)

    deltax_e = n_2(1) - n_1(1);
    deltay_e = n_2(2) - n_1(2);

    L_e = sqrt((deltax_e)^2 + (deltay_e)^2);


    k11 = [A 0 0; 0 12*I/L_e^2 6*I/L_e; 0 6*I/L_e 4*I];

    k12 = [-A 0 0; 0 -12*I/L_e^2 6*I/L_e; 0 -6*I/L_e 2*I];

    k22 = [A 0 0; 0 12*I/L_e^2 -6*I/L_e; 0 -6*I/L_e 4*I];

    k21 = k12';

    Ke = E/L_e*[k11 k12; k21 k22];

    alfa = atan(deltay_e/deltax_e);

    R = [cos(alfa) sin(alfa) 0; -sin(alfa) cos(alfa) 0; 0 0 1];

    Te = [R zeros(3); zeros(3) R];

    K_e = Te'*Ke*Te;  



end

