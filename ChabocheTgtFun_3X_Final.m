function res = ChabocheTgtFun_3X_Final(x,C3,delta_e_pl_stab_vec,delta_sigma_stab_vec,area_stab_vec)

% 2 CICLI A DIVERSO R E DIVERSO DELTA E_PL
% FUNZIONE DI ERRORE DEFINITIVA, OCCHIO AL COEFFICIENTE DI area_stab_vec

% x(1) = C1
% x(2) = g1
% x(3) = C2
% x(4) = g2
% x(5) = sigma0

nCyc = length(delta_e_pl_stab_vec);

res = zeros(2*nCyc,1);

res(1:nCyc) = (2.*(x(1)./x(2)).*(1-exp(-x(2).*delta_e_pl_stab_vec))./(1+exp(-x(2).*delta_e_pl_stab_vec))...
    + 2.*(x(3)./x(4)).*(1-exp(-x(4).*delta_e_pl_stab_vec))./(1+exp(-x(4).*delta_e_pl_stab_vec))...
    + 2.*x(5) - delta_sigma_stab_vec + C3.*delta_e_pl_stab_vec)./(1.*delta_sigma_stab_vec); % Normalizzata a 1

res(nCyc+1:2*nCyc) = (2.*(x(1)./x(2)).*delta_e_pl_stab_vec - 4.*(x(1)./(x(2).^2)).*(1-exp(-x(2).*delta_e_pl_stab_vec))./(1+exp(-x(2).*delta_e_pl_stab_vec))...
    + 2.*(x(3)./x(4)).*delta_e_pl_stab_vec - 4.*(x(3)./(x(4).^2)).*(1-exp(-x(4).*delta_e_pl_stab_vec))./(1+exp(-x(4).*delta_e_pl_stab_vec))...
    + 2.*x(5).*delta_e_pl_stab_vec - area_stab_vec)./(2.*area_stab_vec); % Normalizzata a 1

end

