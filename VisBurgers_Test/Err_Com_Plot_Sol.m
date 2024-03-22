d = 1e-3; % viscosity
m = 400; tf = 1; 
x = linspace(0,1,m+2)';
xi = x(2:end-1);
dx = x(2)-x(1);

CFL = 0.87;
dt = CFL*dx/1;


num_u_331 = compute_sol_visVurgers(d,m,tf,dt,3,3,1,1);
num_u_533 = compute_sol_visVurgers(d,m,tf,dt,5,3,3,3);
true_u = ue(x,tf,d);

% Error
%---------------------------------------------------------%
L1_uerr_331 = dx*sum(abs(true_u-num_u_331));
L1_uerr_533 = dx*sum(abs(true_u-num_u_533));

fprintf('The L_1 error by (3,3,1)=%1.2e \n', L1_uerr_331)
fprintf('The L_1 error by (5,3,3)=%1.2e \n', L1_uerr_533)

L2_uerr_331 = (dx*sum((abs(true_u-num_u_331)).^2))^(1/2);
L2_uerr_533 = (dx*sum((abs(true_u-num_u_533)).^2))^(1/2);

fprintf('The L_2 error by (3,3,1)=%1.2e \n', L2_uerr_331)
fprintf('The L_2 error by (5,3,3)=%1.2e \n', L2_uerr_533)

Linf_uerr_331 = max(abs(true_u-num_u_331));
Linf_uerr_533 = max(abs(true_u-num_u_533));

fprintf('The L inf error by (3,3,1)=%1.2e \n', Linf_uerr_331)
fprintf('The L inf error by (5,3,3)=%1.2e \n', Linf_uerr_533)
%---------------------------------------------------------%

% plot
plot(x,num_u_331,'-r',x,num_u_533,'--b','LineWidth',3)
hold on
plot(x,true_u,':k','LineWidth',3)
legend('(3,3,1)-u','(5,3,3)-u','Exact u','Location','best')
xlabel('x');
ylabel('u', 'Rotation', 0)
set(gca,'FontSize',20)
grid minor
title(sprintf('Solution comparison at t=%d', tf))



%------------------------Function Routine---------------------------------%
function num_u = compute_sol_visVurgers(d,m,tf,dt,s,p,q,scheme_no)
    [A,b,c] = RK_Butcher_Tableau(s,p,q,scheme_no);
    nt = ceil(tf/dt); dt = tf/nt;
    
    x = linspace(0,1,m+2)';
    xi = x(2:end-1);
    h = x(2)-x(1);
    e = ones(m,1);
    
    % 1st and 2nd order spatial derivatives approximations. The first order 
    % spatial derivative is approximated the a third-order upwind biad method
    % as we have steep gradient in the solution, so we want to apply 1st
    % derivative to u^2/2 instead of u. The 2nd derivative is approximated by
    % the second order centered finite difference methods. Look at the book by
    % Hundsdorfer for details.
    
    % M is for 2nd derivative
    M = spdiags([e -2*e e],-1:1,m,m);
    M = M./h^2; M =sparse(M);
    % N is for 1st derivative
    N = spdiags([e -6*e 3*e 2*e],-2:1,m,m); 
    N(1,1) = 2; %Due to the boundary treatment.Also have boundary contributions.
    N = N./(6*h); N =sparse(N);
    
    %---------------------------------------------------------------------%
    s = length(c);
    tn = 0.0;
    u = ue(xi,tn,d); % Inital condition evaluation at tn = 0
    [bc_uxx,bc_fux] = G(xi,tn,d,h);
    
    % Time loop
    for i = 1:nt
        t = tn+dt*c;
        g = zeros(m,s);     % for storing intermediate stages as column vectors
  	    R = zeros(m,s);                      % Evauation of rhs function at g_j
        for j = 1:s
            [bc_uxx,bc_fux] = G(xi,t(j),d,h);
            % Compute the stages
            g_j = u;  
            if j>1
                g_j = g_j + dt*R(:,1:j-1)*(A(j,1:j-1))';
            end
            % Calculating the RHS at (tn+dt*c_j,g_j)
            % R(t,U(t)) = -[NU(t)^2/2 + bc_fux] + d*[MU(t)+bc_uxx] 
            R(:,j) = - ( N*(g_j.^2/2) + bc_fux ) + d*(M*g_j + bc_uxx);
            g(:,j) = g_j;
        end
        tn = tn+dt; % time update
        u = u + dt*R*b'; % solution update
    
        % Plot solution as movie
        % figure(1)
        % clf
        % plot(xi,u,'--b',xi,ue(xi,tn,d),':k','LineWidth',3)
        % drawnow
    end
    
    % Solution at fnal time
    num_u = [ue(x(1),tn,d);u;ue(x(end),tn,d)];
    %true_u = ue(x,tf,d);
    
    % Plot solutions
    % plot(x,num_u,'--b',x,true_u,':k','LineWidth',3)
    % legend('Numerical','Exact')
    
    % Max error
    % uerr = max(abs(true_u-num_u));
end

% Exact Solution
function z = ue(x,t,d)
    r1 = exp(-(20*(x-0.5)+99*t)/(400*d));
    r2 = exp(-(4*(x-0.5)+3*t)/(16*d));
    r3 = exp(-(x-3/8)/(2*d));
    z = 1-0.9*r1./(r1+r2+r3)-0.5*r2./(r1+r2+r3);
end
% Flux function
function z = fue(x,t,d)
    z = 0.5*(ue(x,t,d)).^2;
end
% Exact boundary ondition as vectors
function [bc_uxx,bc_fux] = G(xi,t,d,h)
    bc_uxx = zeros(length(xi),1);
    bc_fux = zeros(length(xi),1);
    % Boundary contribution from 2nd derivative
    bc_uxx(1,1) = g0(t,d)/h^2; bc_uxx(end,1) = g1(t,d)/h^2;
    % Boundary contribution from 1st derivative
    bc_fux(1,1) = -4*fg0(t,d)/(6*h); bc_fux(2,1) = fg0(t,d)/(6*h);
    bc_fux(end,1) = 2*fg1(t,d)/(6*h);
end
% Exact boundary condition left
function z = g0(t,d)
    z = ue(0,t,d);
end
% Exact boundary condition right
function z = g1(t,d)
    z = ue(1,t,d);
end
% Exact boundary condition to approximate (u^2/2)_x: left
function z = fg0(t,d)
    z = fue(0,t,d);
end
% Exact boundary condition to approximate (u^2/2)_x: right
function z = fg1(t,d)
    z = fue(1,t,d);
end
%-------------------------------------------------------------------------%
