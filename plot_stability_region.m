function plot_stability_region(ax,class,varargin)
%==========================================================================
% This code plots the stability region for a given matrix A or (A,b)
% %--Inputs--%
  % ax = axis length
  % Expected class: 'stiff', 'non-stiff'
  % (A,b) or A
% %--Output--%  
  % Plot of stability region.
%==========================================================================
    if strcmp(class,'non-stiff')
        A = varargin{1};
        b = varargin{2};
    elseif strcmp(class,'stiff')
        A = varargin{1};
        b = A(end,:)';
    end
    s = length(b); % number of stages
    e = ones(s,1); % vector of all ones
    I = eye(s); % identity matrix
    zr = linspace(-1,1,200)*ax; zi = zr;
    [ZR,ZI] = meshgrid(zr,zi);
    R = ZR*0;
    for jr = 1:length(ZR)
        for ji = 1:length(ZI)
            z = zr(jr)+1i*zi(ji);
            R(ji,jr) = 1+z*b'*((I-z*A)\e);
        end
    end
    contourf(zr,zi,1-abs(R),[0 0])
    hold on
    xline(0)
    hold on
    yline(0)
    % title('Stability region')
    axis equal
    set(gca,'FontSize',30)
    % to have minimum white space
     ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end