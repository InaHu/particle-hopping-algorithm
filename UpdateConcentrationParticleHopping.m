function [a1, r1, j] = UpdateConcentrationParticleHopping(N) 
    N.p0 = N.a0 + N.r0; % Sum of concentrations
    
    C_a = N.Scaling_Drift/(2*N.eps_a*N.Scaling_Diffusion);
    C_r = N.Scaling_Drift/(2*N.eps_r*N.Scaling_Diffusion);
    
    a1 = zeros(length(N.a0),1); % New concentrations
    r1 = zeros(length(N.r0),1);
    
    G = - N.a0(2:end-1).*(1-N.p0(1:end-2)).*exp(-C_a*(N.V_a(2:end-1)-N.V_a(1:end-2)));          % Probability that a ANT-particle jumps away to the left
    G = G + N.a0(1:end-2).*(1-N.p0(2:end-1)).*exp(-C_a*(N.V_a(1:end-2)-N.V_a(2:end-1)));
    G = G - N.a0(2:end-1).*(1-N.p0(3:end)).*exp(-C_a*(N.V_a(2:end-1)-N.V_a(3:end)));
    G = G + N.a0(3:end).*(1-N.p0(2:end-1)).*exp(-C_a*(N.V_a(3:end)-N.V_a(2:end-1)));
    
    a1(2:end-1) = N.a0(2:end-1) + N.Scaling_Diffusion*(N.eps_a*N.dt/N.dx^2)*G;
    
    G = - N.r0(2:end-1).*(1-N.p0(1:end-2)).*exp(-C_r*(N.V_r(2:end-1)-N.V_r(1:end-2)));          % Probability that a RET particle jumps away to the left
    G = G + N.r0(1:end-2).*(1-N.p0(2:end-1)).*exp(-C_r*(N.V_r(1:end-2)-N.V_r(2:end-1)));
    G = G - N.r0(2:end-1).*(1-N.p0(3:end)).*exp(-C_r*(N.V_r(2:end-1)-N.V_r(3:end)));
    G = G + N.r0(3:end).*(1-N.p0(2:end-1)).*exp(-C_r*(N.V_r(3:end)-N.V_r(2:end-1)));
    
    r1(2:end-1) = N.r0(2:end-1) + N.Scaling_Diffusion*(N.eps_r*N.dt/N.dx^2)*G;
    
    C_a = N.Scaling_Drift/(2*N.eps_a*N.Scaling_Diffusion) * N.typicalLength/N.typicalTime ;
    C_r = N.Scaling_Drift/(2*N.eps_r*N.Scaling_Diffusion) * N.typicalLength/N.typicalTime;
    
    % Boundary Flow at x = 0
    a1(1) = N.a0(1) + N.eps_a*N.dt/(N.dx)^2*(-N.a0(1)*(1-N.p0(2))*exp(-C_a*(N.V_a(1)-N.V_a(2))) ... 
        + N.a0(2)*(1-N.p0(1))*exp(-C_a*(N.V_a(2)-N.V_a(1))) ... 
        + C_a*N.dx*N.Lambda_som*N.alpha_a*(1-N.p0(1)));
    
    r1(1) = N.r0(1) + N.eps_r*N.dt/(N.dx)^2*(-N.r0(1)*(1-N.p0(2))*exp(-C_r*(N.V_r(1)-N.V_r(2))) ...
        + N.r0(2)*(1-N.p0(1))*exp(-C_r*(N.V_r(2)-N.V_r(1))) ...
        - C_r*N.dx*N.beta_r*N.r0(1)*(N.Lambda_som_max - N.Lambda_som));
    
    % Boundary Flow at x = 1 
    a1(end) = N.a0(end) ...
        + N.eps_a*N.dt/(N.dx)^2*(-N.a0(end)*(1-N.p0(end-1))*exp(-C_a*(N.V_a(end)-N.V_a(end-1))) ...
        + N.a0(end-1)*(1-N.p0(end))*exp(-C_a*(N.V_a(end-1)-N.V_a(end))) ... 
        - C_a*N.dx*N.beta_a*N.a0(end)*(N.Lambda_tip_max - N.Lambda_tip));
    
    r1(end) = N.r0(end) ...
        + N.eps_r*N.dt/(N.dx)^2*(-N.r0(end)*(1-N.p0(end-1))*exp(-C_r*(N.V_r(end)- N.V_r(end-1))) ... 
        + N.r0(end-1)*(1-N.p0(end))*exp(-C_r*(N.V_r(end-1)-N.V_r(end))) ... 
        + C_r*N.dx*N.alpha_r*N.Lambda_tip*(1-N.p0(end)));
    
    j = N.a0 + N.r0 - a1 - r1;
end