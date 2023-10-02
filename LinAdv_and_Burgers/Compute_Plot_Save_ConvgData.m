%=========================================================================%
% This code computes the data for convergence plots for the linear        %
% advection equation and the Burgers' equation. Change the parameter      %
% 'eq_ind' for different test cases below.                                %
%=========================================================================%

eq_ind = 2; % change here for different test case
% Different test problems
EQN = {'LinAdv','Burgers'};
% Different initial conditions
IC = [1,1];
% Final time
TF = [0.7,0.8];

eqn=EQN{eq_ind}; TC=IC(eq_ind); tf=TF(eq_ind);

switch eqn
    case 'LinAdv'
        % Spatial order and grid points
        sp_or = [1,1,1,1,1,1,1,1,1];
        max_a = 2; % max of the inital profile for the Burgers equation
        M = [ceil(10.^linspace(1.2,3,10));ceil(10.^linspace(1.2,3,10));ceil(10.^linspace(1.2,3,10));
            ceil(10.^linspace(1.2,2.5,10));ceil(10.^linspace(1.2,2.5,10));ceil(10.^linspace(1.2,2.5,10));
            ceil(10.^linspace(.8,2.2,10));ceil(10.^linspace(.8,2.2,10));ceil(10.^linspace(.8,2.2,10))];
        %M = ceil(10.^linspace(1.2,3,10));
        MU = [0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9]; % CFL numbers
    case 'Burgers'
        % Spatial order and grid points
        sp_or = [1,1,1,1,1,1,1,1,1];
        max_a = 2; % max of the inital profile for the Burgers equation
        M = [ceil(10.^linspace(1.2,3,10));ceil(10.^linspace(1.2,3,10));ceil(10.^linspace(1.2,3,10));
            ceil(10.^linspace(1.2,2.5,10));ceil(10.^linspace(1.2,2.5,10));ceil(10.^linspace(1.2,2.5,10));
            ceil(10.^linspace(1.2,2.5,10));ceil(10.^linspace(1.2,2.5,10));ceil(10.^linspace(1.2,2.5,10))];
        %M = ceil(10.^linspace(1.2,3,10));
        MU = [0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9]; % CFL numbers
end

% ERK-(s,p,q) scheme
S=[3,4,5,4,6,7,7,8,9]; 
P=[3,3,3,4,4,4,5,5,5]; 
Q=[1,2,3,1,3,4,1,4,5]; 
SchNo=[1,2,3,4,5,6,7,8,9];

% Compute

U_Err = zeros(9,length(M)); dU_Err = zeros(9,length(M));
for i = 1:length(S)
    mu = MU(i);
    H = 1./M; dts = mu*H/max_a;
    s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); spatial_order = sp_or(i);    
    parfor j =1:length(M(i,:))
        switch eqn
            case 'LinAdv'
                [uerr,duerr] = LinAdvTestProbZeroSpError(TC,M(i,j),s,p,q,scheme_no,mu,tf);
                U_Err(i,j) = uerr; dU_Err(i,j) = duerr; 
            case 'Burgers'
                [uerr,duerr] = BurgersTestProbZeroSpError(TC,M(i,j),s,p,q,scheme_no,mu,tf);
                U_Err(i,j) = uerr; dU_Err(i,j) = duerr;              
        end
    end   
end

% saving the data
filename = sprintf('%s_ConvgData.mat',eqn);
save(filename, 'dts', 'U_Err', 'dU_Err');

%plot
C = {'b','r','g'}; Cref = {[0.5,0.5,0.5]};
linS = {'-',':'}; Mar = {'o','s','+'}; ms = 8; fs = 12;

% for legends
switch eqn 
    case 'LinAdv'
        Cof3 = [1e-1,9e-2,2e+0]; Sl3 = [1,2,3]; lim3 = [1e-11,5e-3];
        st31 = 7; end31 = 10; st32 = 5; end32 = 8; st33 = 7; end33 = 10;
        Cof4 = [5e-3,6e-2,3e+0]; Sl4 = [2,3,4]; lim4 = [5e-15,1e-3];
        st41 = 7; end41 = 10; st42 = 6; end42 = 10; st43 = 4; end43 = 7;
        Cof5 = [1e-4,3e-2,5e-2]; Sl5 = [2,4,5]; lim5 = [5e-16,1e-5];
        st51 = 7; end51 = 10; st52 = 3; end52 = 6; st53 = 1; end53 = 4;
    case 'Burgers'
        Cof3 = [8e-3,5e-3,5e-1]; Sl3 = [1,2,3]; lim3 = [8e-12,1e-2];
        st31 = 7; end31 = 10; st32 = 7; end32 = 10;  st33 = 7; end33 = 10;
        Cof4 = [1e-3,4e-3,5e-3]; Sl4 = [1,2,3]; lim4 = [1e-14,1e-4]; 
        st41 = 7; end41 = 10; st42 = 7; end42 = 10; st43 = 7; end43 = 10;
        Cof5 = [8e-6,8e-4,7e-5]; Sl5 = [1,2,3];lim5 = [5e-15,5e-5];
        st51 = 5; end51 = 10; st52 = 5; end52 = 8; st53 = 5; end53 = 8;
end


%-------------------------------%
% save_data = 1;
% if save_data 
%     foldername_Figures = sprintf('Figures/');
%     if exist(foldername_Figures,'dir') == 0
%         mkdir(foldername_Figures);
%     end
% end
%-------------------------------%
figure(1)
set(0,'DefaultLineLineWidth',3);
set(gcf,'position',[50 50 400 550])
axis_position = [.18 .18 .80 .80];
set(gca,'position',axis_position)
set(gcf,'paperpositionmode','auto')
set(gcf,'PaperPosition',[0 0 400/100 550/100])
for i = 1:3
    loglog(dts(i,:),U_Err(i,:),'linestyle',linS{1},'color',C{i},'marker',Mar{i},'MarkerSize',ms)
    hold on
    loglog(dts(i,:),dU_Err(i,:),'linestyle',linS{2},'color',C{i},'marker',Mar{i},'MarkerSize',ms)
    hold on
end
hold on
loglog(dts(1,st31:end31),Cof3(1)*dts(1,st31:end31).^Sl3(1),'--','color',Cref{1})
loglog(dts(1,st32:end32),Cof3(2)*dts(1,st32:end32).^Sl3(2),'-.','color',Cref{1})
loglog(dts(1,st33:end33),Cof3(3)*dts(1,st33:end33).^Sl3(3),'-','color',Cref{1})
legend('','Location', 'northoutside')
switch eqn
    case 'LinAdv'
        legend({'(3,3,1):u','(3,3,1):u_x','(4,3,2):u','(4,3,2):u_x',...
        '(5,3,3):u','(5,3,3):u_x',sprintf('Slope %d',Sl3(1)),...
        sprintf('Slope %d',Sl3(2)),sprintf('Slope %d',Sl3(3))},'NumColumns',3)
    case 'Burgers'
        legend({'(3,3,1):u','(3,3,1):u_x','(4,3,2):u','(4,3,2):u_x',...
        '(5,3,3):u','(5,3,3):u_x',sprintf('Slope %d',Sl3(1)),...
        sprintf('Slope %d',Sl3(2)),sprintf('Slope %d',Sl3(3))},'NumColumns',3)
end
xlim([dts(1,end),dts(1,1)])
ylim([lim3(1),lim3(2)])
xlabel('\Delta t');
ylabel('Error');
set(gca,'FontSize',fs)

% Save in eps
% figure_name = sprintf('Figures/%s_ConvPlotOrder_%d',eqn,P(1));
% print(gcf,figure_name,'-depsc','-r100','-vector')
%-------------------------------%
figure(2)
set(0,'DefaultLineLineWidth',3);
set(gcf,'position',[50 50 400 550])
axis_position = [.18 .18 .80 .80];
set(gca,'position',axis_position)
set(gcf,'paperpositionmode','auto')
set(gcf,'PaperPosition',[0 0 400/100 550/100])
for i = 4:6
    loglog(dts(i,:),U_Err(i,:),'linestyle',linS{1},'color',C{i-3},'marker',Mar{i-3},'MarkerSize',ms)
    hold on
    loglog(dts(i,:),dU_Err(i,:),'linestyle',linS{2},'color',C{i-3},'marker',Mar{i-3},'MarkerSize',ms)
    hold on
end
hold on
loglog(dts(4,st41:end41),Cof4(1)*dts(4,st41:end41).^Sl4(1),'--','color',Cref{1})
loglog(dts(4,st42:end42),Cof4(2)*dts(4,st42:end42).^Sl4(2),'-.','color',Cref{1})
loglog(dts(4,st43:end43),Cof4(3)*dts(4,st43:end43).^Sl4(3),'-','color',Cref{1})

legend('','Location', 'northoutside')
switch eqn
    case 'LinAdv'
        legend({'(4,4,1):u','(4,4,1):u_x','(6,4,3):u','(6,4,3):u_x',...
        '(7,4,4):u','(7,4,4):u_x',sprintf('Slope %d',Sl4(1)),...
        sprintf('Slope %d',Sl4(2)),sprintf('Slope %d',Sl4(3))},'NumColumns',3)
    case 'Burgers'
        legend({'(4,4,1):u','(4,4,1):u_x','(6,4,3):u','(6,4,3):u_x',...
        '(7,4,4):u','(7,4,4):u_x',sprintf('Slope %d',Sl4(1)),...
        sprintf('Slope %d',Sl4(2)),sprintf('Slope %d',Sl4(3))},'NumColumns',3)
end
xlim([dts(4,end),dts(4,1)])
ylim([lim4(1),lim4(2)])
xlabel('\Delta t');
ylabel('Error');
set(gca,'FontSize',fs)

% Save in eps
% figure_name = sprintf('Figures/%s_ConvPlotOrder_%d',eqn,P(4));
% print(gcf,figure_name,'-depsc','-r100','-vector')
%-------------------------------%
figure(3)
set(0,'DefaultLineLineWidth',3);
set(gcf,'position',[50 50 400 550])
axis_position = [.18 .18 .80 .80];
set(gca,'position',axis_position)
set(gcf,'paperpositionmode','auto')
set(gcf,'PaperPosition',[0 0 400/100 550/100])
for i = 7:9
    loglog(dts(i,:),U_Err(i,:),'linestyle',linS{1},'color',C{i-6},'marker',Mar{i-6},'MarkerSize',ms)
    hold on
    loglog(dts(i,:),dU_Err(i,:),'linestyle',linS{2},'color',C{i-6},'marker',Mar{i-6},'MarkerSize',ms)
    hold on
end
hold on
loglog(dts(7,st51:end51),Cof5(1)*dts(7,st51:end51).^Sl5(1),'--','color',Cref{1})
loglog(dts(7,st52:end52),Cof5(2)*dts(7,st52:end52).^Sl5(2),'-.','color',Cref{1})
loglog(dts(7,st53:end53),Cof5(3)*dts(7,st53:end53).^Sl5(3),'-','color',Cref{1})
legend('','Location', 'northoutside')
switch eqn
    case 'LinAdv'
        legend({'(7,5,1):u','(7,5,1):u_x','(8,5,4):u','(8,5,4):u_x',...
        '(9,5,5):u','(9,5,5):u_x',sprintf('Slope %.1f',Sl5(1)),...
        sprintf('Slope %d',Sl5(2)),sprintf('Slope %d',Sl5(3))},'NumColumns',3)
    case 'Burgers'
        legend({'(7,5,1):u','(7,5,1):u_x','(8,5,4):u','(8,5,4):u_x',...
        '(9,5,5):u','(9,5,5):u_x',sprintf('Slope %d',Sl5(1)),...
        sprintf('Slope %d',Sl5(2)),sprintf('Slope %d',Sl5(3))},'NumColumns',3)
end
xlim([dts(7,end),dts(7,1)])
ylim([lim5(1),lim5(2)])
xlabel('\Delta t');
ylabel('Error');
set(gca,'FontSize',fs)

% Save in eps
% figure_name = sprintf('Figures/%s_ConvPlotOrder_%d',eqn,P(7));
% print(gcf,figure_name,'-depsc','-r100','-vector')
%-------------------------------%


