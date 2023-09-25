function [A,b,c,exitflag] = FindOptSchemes(s,p,q,no_it)
%=========================================================================%
%This code finds explicit Runge-Kutta methods (outputs: A,b) for the given%
% parameters (inputs: s, p, q, no_it).                                    %
%=========================================================================%     
    % Define the objective function
    obj_fun = @(x) 0; 
    %obj_fun = @(x) objective_fun(x, s, p);
    
    % Define the nonlinear constraints
    nonlconfun = @(x) myconstraints(x, s, p, q);
    
    % Define the lower and upper bounds for x
    max_mag = 7; % This can be set to inf if just testing whether methods exist
    lb = [-max_mag*ones(s*(s-1)/2, 1); -max_mag*ones(s, 1)];
    ub = [max_mag*ones(s*(s-1)/2, 1); max_mag*ones(s, 1)];
    
    % Solve the optimization problem
    opts = optimoptions(@fmincon,'MaxFunEvals',20000,'TolCon',1.e-14,...
        'TolFun',1.e-14,'TolX',1.e-14,'MaxIter',20000,'Algorithm','sqp',...
        'Display','notify');
    exitflag = -1;
    iter = 0;
    
    while (exitflag ~= 1) && (iter<no_it)
        % Define the initial guess
        A0 = tril(rand(s),-1); b0 = rand(s,1); x0 = pack_rk(A0,b0);
    
        [x,fval,exitflag,output,lambda] = fmincon(obj_fun,x0,[],[],[],[],...
            lb,ub, nonlconfun,opts);
        iter = iter+1
    end
    
    % Extract A and b from x
    [A,b] = unpack_rk(x,s);
    c = sum(A,2);
    [Cineq, Ceq] = myconstraints(x, s, p, q);
    
    fprintf('Exitflag = %d \n',exitflag)
    
    % save the method
    % if exitflag == 1
    %     MthdName = sprintf('ExMthd_s%dp%dq%d.mat',s,p,q);
    %     save(MthdName,'A','b')
    % end
end
%=========================================================================%
function [Cineq, Ceq] = myconstraints(x, s, p, q)
    % No inequality constraints
    Cineq = [];
    % Extract A and b from x
    [A,b] = unpack_rk(x,s);

    % Compute the values of e, c, and tau(k)
    e = ones(s,1); c = A*e; C = diag(c); Tau = residual(A,c,s,q);
    % WSO conditons
    Ceq = zeros();
    for i = 1:q-1
        for j = 1:s
            Ceq((i-1)*s+j,1) = b'*A^(j-1)*Tau(:,i);
        end
    end
    % 1st order conditons
    if p>=1
        Ceq((q-1)*s+1,1) = b'*e-1;
    end
    % 2nd order conditons
    if p>=2
        Ceq((q-1)*s+2,1) = b'*c-1/2;
    end
    % 3rd order conditons
    if p>=3
        Ceq((q-1)*s+3,1) = b'*A*c-1/6;
        Ceq((q-1)*s+4,1) = b'*c.^2-1/3;
    end
    % 4th order conditons
    if p>=4  
        Ceq((q-1)*s+5,1) = b'*C*A*c-1/8;
        Ceq((q-1)*s+6,1) = b'*A^2*c-1/24;
        Ceq((q-1)*s+7,1) = b'*A*c.^2-1/12;
        Ceq((q-1)*s+8,1) = b'*c.^3-1/4;
    end
    if p>=5
        %5th order conditions
        Ceq((q-1)*s+9,1) = b'*c.^4-1/5;
        Ceq((q-1)*s+10,1) = b'*((A*c).*c.^2)-1/10;
        Ceq((q-1)*s+11,1) = b'*((A*c.^2).*c)-1/15;
        Ceq((q-1)*s+12,1) = b'*(A*c).^2-1/20;
        Ceq((q-1)*s+13,1) = b'*((A^2*c).*c)-1/30;
        Ceq((q-1)*s+14,1) = b'*A*C*A*c-1/40;
        Ceq((q-1)*s+15,1) = b'*A*c.^3-1/20;
        Ceq((q-1)*s+16,1) = b'*A^2*c.^2-1/60;
        Ceq((q-1)*s+17,1) = b'*A^3*c-1/120;
    end
    if p>=6
        %6th order conditions
        Ceq((q-1)*s+18,1) = c'.^5*b-1/6;
        Ceq((q-1)*s+19,1) = b'*diag(c).^3*A*c-1/12;
        Ceq((q-1)*s+20,1) = b'*diag(c)*(A*c).^2-1/24;
        Ceq((q-1)*s+21,1) = b'*diag(c).^2*A*c.^2-1/18;
        Ceq((q-1)*s+22,1) = b'*((A*c.^2).*(A*c))-1/36;
        Ceq((q-1)*s+23,1) = b'*diag(c)*A*c.^3-1/24;
        Ceq((q-1)*s+24,1) = b'*A*c.^4-1/30;
        Ceq((q-1)*s+25,1) = b'*diag(c).^2*A^2*c-1/36;
        Ceq((q-1)*s+26,1) = b'*((A^2*c).*(A*c))-1/72;
        Ceq((q-1)*s+27,1) = b'*diag(c)*A*diag(c)*A*c-1/48;
        Ceq((q-1)*s+28,1) = b'*A*diag(c).^2*A*c-1/60;
        Ceq((q-1)*s+29,1) = b'*A*(A*c).^2-1/120;
        Ceq((q-1)*s+30,1) = b'*diag(c)*A^2*c.^2-1/72;
        Ceq((q-1)*s+31,1) = b'*A*diag(c)*A*c.^2-1/90;
        Ceq((q-1)*s+32,1) = b'*A^2*c.^3-1/120;
        Ceq((q-1)*s+33,1) = b'*diag(c)*A^3*c-1/144;
        Ceq((q-1)*s+34,1) = b'*A*diag(c)*A^2*c-1/180;
        Ceq((q-1)*s+35,1) = b'*A^2*diag(c)*A*c-1/240;
        Ceq((q-1)*s+36,1) = b'*A^3*c.^2-1/360;
        Ceq((q-1)*s+37,1) = b'*A^4*c-1/720;           
    end
    Ceq(end+1,1) = A(q+1,q); % Butcher tableau structural restriction
end

%-------------------------------------------------------------------------%
% Stage order residuals.
%-------------------------------------------------------------------------%
function tau = residual(A,c,s,q)  
    tau = zeros(s,q-1);
     for i = 1:q-1
         tau(:,i) = c.^(i+1) - (i+1)*A*c.^(i); 
     end
end

%-------------------------------------------------------------------------%
% This function returns A and b from the long vector x.
%-------------------------------------------------------------------------%
function [A,b] = unpack_rk(x,s) 
    A = zeros(s,s); 
     for i = 2:s
         A(i,1:i-1) = x(((i-1)*(i-2)/2+1):(i*(i-1)/2));
     end 
     b = x(s*(s-1)/2+1:s*(s-1)/2+s,1);
end

%-------------------------------------------------------------------------%
% This function takes A,b and reshapes them into a long vector.
%-------------------------------------------------------------------------%
function x = pack_rk(A,b)
    s = size(A,1);
    x = zeros((s*(s-1)/2+s),1);
    for i = 2:s
        x( ((i-1)*(i-2)/2+1):(i*(i-1)/2),1 )=A(i,1:i-1)';
    end
    x(s*(s-1)/2+1:s*(s-1)/2+s,1) = b;
end
%-------------------------------------------------------------------------%
% The objective function.
%-------------------------------------------------------------------------%
function y = objective_fun(x, s, p)
    % Extract A and b from x
    [A,b] = unpack_rk(x,s);
    % Compute the values of e, c, and tau(k)
    e = ones(s,1); c = A*e; C = diag(c); 
    Ceq = zeros();
    if p == 4
        %5th order conditions
        Ceq(1) = b'*c.^4-1/5;
        Ceq(2) = b'*((A*c).*c.^2)-1/10;
        Ceq(3) = b'*((A*c.^2).*c)-1/15;
        Ceq(4) = b'*(A*c).^2-1/20;
        Ceq(5) = b'*((A^2*c).*c)-1/30;
        Ceq(6) = b'*A*C*A*c-1/40;
        Ceq(7) = b'*A*c.^3-1/20;
        Ceq(8) = b'*A^2*c.^2-1/60;
        Ceq(9) = b'*A^3*c-1/120;
    end
    y = norm(Ceq);
end
%==================================End====================================%



