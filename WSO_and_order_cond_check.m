function [OC_err,WSO_err] = WSO_and_order_cond_check(s,p,q,class,varargin)
%==========================================================================
% This code checks whether a pair (A,b) satisfies pth (<=6) order conditions
% and the weak stage order q (>=2) conditions.
% %--Inputs--%
  % s = Any number of stages
  % Expected class: 'stiff', 'non-stiff'
  % (A,b) or A
% %--Output--%  
  % error in order conditions and WSO conditions.
%==========================================================================
    if strcmp(class,'non-stiff')
        A = varargin{1};
        b = varargin{2};
    elseif strcmp(class,'stiff')
        A = varargin{1};
        b = A(end,:)';
    end
    e= ones(s,1);
    c = sum(A,2);
    C = diag(c); % diagonal matrix of the vector c
    res = zeros(s,q-1); % storing residual vectors r2, r3, and r4 in the columns of res
    for i = 1:q-1
        res(:,i) = A*c.^(i)-(1/(i+1))*c.^(i+1);
    end

    if p>=1
        % 1st order condition
        OC_err(1) = b'*e-1;
    end
    if p>=2
        % 2nd order condition
        OC_err(2) = b'*c-1/2;
    end
    if p>=3
        % 3rd order conditions
        OC_err(3) = b'*A*c-1/6;
        OC_err(4) = b'*c.^2-1/3;
    end
    if p>=4  
        % 4th order conditons
        OC_err(5) = b'*A^2*c-1/24;
        OC_err(6) =b'*A*c.^2-1/12;
        OC_err(7) = b'*C*A*c-1/8;
        OC_err(8) = b'*c.^3-1/4;
    end
    if p>=5
        %5th order conditions
        OC_err(9) = b'*c.^4-1/5;
        OC_err(10) = b'*((A*c).*c.^2)-1/10;
        OC_err(11) = b'*((A*c.^2).*c)-1/15;
        OC_err(12) = b'*(A*c).^2-1/20;
        OC_err(13) = b'*((A^2*c).*c)-1/30;
        OC_err(14) = b'*A*C*A*c-1/40;
        OC_err(15) = b'*A*c.^3-1/20;
        OC_err(16) = b'*A^2*c.^2-1/60;
        OC_err(17) = b'*A^3*c-1/120;
    end
    if p>=6
        %6th order conditions
        OC_err(18)=c'.^5*b-1/6;
        OC_err(19)=b'*diag(c).^3*A*c-1/12;
        OC_err(20)=b'*diag(c)*(A*c).^2-1/24;
        OC_err(21)=b'*diag(c).^2*A*c.^2-1/18;
        OC_err(22)=b'*((A*c.^2).*(A*c))-1/36;
        OC_err(23)=b'*diag(c)*A*c.^3-1/24;
        OC_err(24)=b'*A*c.^4-1/30;
        OC_err(25)=b'*diag(c).^2*A^2*c-1/36;
        OC_err(26)=b'*((A^2*c).*(A*c))-1/72;
        OC_err(27)=b'*diag(c)*A*diag(c)*A*c-1/48;
        OC_err(28)=b'*A*diag(c).^2*A*c-1/60;
        OC_err(29)=b'*A*(A*c).^2-1/120;
        OC_err(30)=b'*diag(c)*A^2*c.^2-1/72;
        OC_err(31)=b'*A*diag(c)*A*c.^2-1/90;
        OC_err(32)=b'*A^2*c.^3-1/120;
        OC_err(33)=b'*diag(c)*A^3*c-1/144;
        OC_err(34)=b'*A*diag(c)*A^2*c-1/180;
        OC_err(35)=b'*A^2*diag(c)*A*c-1/240;
        OC_err(36)=b'*A^3*c.^2-1/360;
        OC_err(37)=b'*A^4*c-1/720;           
    end

    WSO_err = zeros(q-1,s); % WSO 4 condition errors
    for i = 1:q-1
        for j = 1:s
            WSO_err(i,j) = b'*A^(j-1)*res(:,i);
        end
    end
end