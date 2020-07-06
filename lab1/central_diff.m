function du = central_diff(u, direction,g)
    M = size(u,1);
    [du_p,rho_p] = dir_diff(u, direction(1), direction(2), g);
    [du_m,rho_m] = dir_diff(u,-direction(1), -direction(2), g);
    du = 2./(direction(1)^2+direction(2)^2).*(M+1)^2./(rho_p + rho_m).*(du_p./rho_p + du_m./rho_m);
end