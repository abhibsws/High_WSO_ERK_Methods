function [A,b,c] = RK_Butcher_Tableau(s,p,q,scheme_no)
%=========================================================================%
% This code gives different RK methods depending on (s,p,q). For tuple 
% (s,p,q), there are a few good schemes and scheme_no determines them
% uniquely. 
% ---Inputs---
    % s = number of stages of RK scheme.
    % p = order of RK scheme.
    % q = weak stage order of RK scheme.
% ---Outputs---
    % Butcher tanleau (A,b,c)
%=========================================================================%
    if s==5 && p ==3 && q==3 && scheme_no==1
        %ExHighWSO-(5,3,3): principal error norm = 0.07. -- our method
        A = [0,0,0,0,0;
            (3/11),0,0,0,0;
            (285645/493487),(103950/493487),0,0,0;
            (3075805/5314896),(1353275/5314896),0,0,0;
            (196687/177710),(-129383023/426077496),(48013/42120),(-2268/2405),0];
        b = [(5626/4725),(-25289/13608),(569297/340200),(324/175),(-13/7)];
        c = [0,(3/11),(15/19),(5/6),1]';
        %ERK313: principal error norm = . -- by Skvortsov
    elseif s==5 && p ==3 && q==3 && scheme_no==2
        A=[     0     0    0   0   0;
              1/3     0    0   0   0;
              2/3     0    0   0   0;
                1     0    0   0   0;
           -11/12   3/2 -3/4 1/6   0];
        b=[1/4,-3,15/4,-1,1];
        c = sum(A,2);
    elseif s==12 && p ==3 && q==3 && scheme_no==3
    % Parallel iterative method: The vector c is the 4-point Gauss quadrature
        c_gauss = [ 0.069431844202974   0.330009478207572   0.669990521792428   0.930568155797026];
        [A, b, c] = pirk(c_gauss);
    end
end
%=========================================================================%








