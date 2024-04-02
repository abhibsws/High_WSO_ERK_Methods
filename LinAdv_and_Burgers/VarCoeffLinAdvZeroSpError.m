function [uerr,duerr] = VarCoeffLinAdvZeroSpError(TC,m,s,p,q,scheme_no,mu,tf) 
    % Spatial approximation by finite difference 
    x = linspace(0,1,m+1)';
    xi = x(2:m+1);
    h = x(2)-x(1);
    e = ones(m,1);
    % N is for u_x
    N = spdiags([e -1*e 0*e],-1:1,m,m);
    N = N./h; N =sparse(N);
    
    % Time stepping 
    s_max = max(abs(vel(x)));
    dt = mu*h/s_max;
    nt = ceil(tf/dt); dt = tf/nt;
    
    % Explicit RK scheme 
    [A,b,c] = RK_Butcher_Tableau(s,p,q,scheme_no);
    
    %-------------------------------------------------------------------------%
    s = length(c);
    tn = 0.0; 
    u = ue(TC,xi,tn); % Inital condition evaluation at tn = 0
    avec = vel(xi); % spatially varying velocity profile
    % Time loop
    for i = 1:nt
        t = tn+dt*c;
        g = zeros(m,s);     % for storing intermediate stages as column vectors
        R = zeros(m,s);                      % Evauation of rhs function at g_j
        for j = 1:s
            g_j = u; 
            if j>1
                g_j= g_j + dt*R(:,1:j-1)*(A(j,1:j-1))';
            end
            bc_ux = G(TC,xi,t(j),h); f_vec = F(TC,xi,t(j));
            % Calculating the right hand side: avec.*(N*U(t) + G(t) + F(t))
            R(:,j) = avec.*(N*g_j + bc_ux) + f_vec;
            g(:,j) = g_j;
        end
        u = u + dt*(R*b');
        tn = tn+dt;
        if i==nt
            ubnd = ue(TC,x(1),tn); num_uvec = [ubnd;u];
        end
    end
    true_uvec = ue(TC,x,tf);
    uerr=max(abs(num_uvec-true_uvec));

    num_duvec = duapprox(h,num_uvec);
    true_duvec = ue_x(TC,x,tf);
    duerr = max(abs(num_duvec-true_duvec));
end
%-------------------------------------------------------------------------%
% Exact Solution
function z = ue(TC,x,t)
    switch TC
        case 1, z = (1+x)/(1+t);
    end
end

%First order time derivative of exact solution
function z = ue_t(TC,x,t)
    switch TC
        case 1, z = -(1+x)/(1+t)^2;
    end
end

%First order spatial derivative of exact solution
function z = ue_x(TC,x,t)
    switch TC
        case 1, z = (0*x+1)/(1+t);
    end
end
function a = vel(x)
    % a = cos(x+0.1);  
    a = (sin(x)+2)/3;
end
% Exact forcing term as vectors
function y = F(TC,x,t)
    y = ue_t(TC,x,t)+vel(x).*ue_x(TC,x,t);
end
% Exact boundary ondition as vectors
function bc_ux = G(TC,xi,t,h)
    bc_ux = zeros(length(xi),1);
    bc_ux(1,1) = 1*g0(TC,t)/h;
end
% Exact boundary condition left (only need this Dirchlet bc as the information
% is moving to the right  with velocity 1).
function z = g0(TC,t)
    z = ue(TC,0,t);
end
% Different order derivative approximation from function values
function du = duapprox(dx,u)
% compute u_x from u: 6th order approximation
du = [sum([-147;360;-450;400;-225;72;-10].*u(1:7))./(60*dx);...
      sum([-10;-77;150;-100;50;-15;2].*u(1:7))./(60*dx);...
      sum([2;-24;-35;80;-30;8;-1].*u(1:7))./(60*dx);...
      (-1*u(1:end-6)+9*u(2:end-5)-45*u(3:end-4)+45*u(5:end-2)-9*u(6:end-1)+1*u(7:end))./(60*dx);...
      sum([1;-8;30;-80;35;24;-2].*u(end-6:end))./(60*dx);...
      sum([-2;15;-50;100;-150;77;10].*u(end-6:end))./(60*dx);...
      sum([10;-72;225;-400;450;-360;147].*u(end-6:end))./(60*dx)];
end
%-------------------------------------------------------------------------%
