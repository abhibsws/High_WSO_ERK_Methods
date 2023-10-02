%=========================================================================%
% This code computes the data for convergence plots for the linear        %
% advection equation and the Burgers' equation. Change the parameter      %
% 'eq_ind' for different test cases below.                                %
%=========================================================================%

eq_ind = 1; % change here for different test case
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
        sp_or = [1,1,1];
        max_a = 2; % max of the inital profile for the Burgers equation
        M = ceil(10.^linspace(1.2,3,10));
        MU = [0.9,0.9,0.9]; % CFL numbers
    case 'Burgers'
        % Spatial order and grid points
        sp_or = [1,1,1];
        max_a = 2; % max of the inital profile for the Burgers equation
        M = ceil(10.^linspace(1.2,3,10));
        MU = [0.9,0.9,0.9]; % CFL numbers
end

% ERK-(s,p,q) scheme
S = [5,5,12]; 
P = [3,3,3]; 
Q = [3,3,3]; 
SchNo=[1,2,3];

% Compute
U_Err = zeros(length(S),length(M)); dU_Err = zeros(length(S),length(M));

for i = 1:length(S)
    mu = MU(i);
    H = 1./M; dts = mu*H/max_a;
    s=S(i); p=P(i); q=Q(i); scheme_no=SchNo(i); spatial_order = sp_or(i);    
    parfor j =1:length(M)
        switch eqn
            case 'LinAdv'
                [uerr,duerr] = LinAdvTestProbZeroSpError(TC,M(j),s,p,q,scheme_no,mu,tf);
                U_Err(i,j) = uerr; dU_Err(i,j) = duerr; 
            case 'Burgers'
                [uerr,duerr] = BurgersTestProbZeroSpError(TC,M(j),s,p,q,scheme_no,mu,tf);
                U_Err(i,j) = uerr; dU_Err(i,j) = duerr;              
        end
    end   
end

% saving the data
filename = sprintf('%s_ConvgData_WSO3.mat',eqn);
save(filename, 'dts', 'U_Err', 'dU_Err');

%plot
C = {'b','r','g'}; Cref = {[0.5,0.5,0.5]};
linS = {'-',':'}; Mar = {'o','s','+'}; ms = 8; fs = 14;

% for legends
switch eqn 
    case 'LinAdv'
        Cof3 = [1e-1,9e-2,2e+0]; Sl3 = [3]; lim3 = [1e-11,5e-5];
        st31 = 5; end31 = 10; 
  case 'Burgers'
        Cof3 = [8e-3,5e-3,5e-1]; Sl3 = [3]; lim3 = [8e-12,1e-5];
        st31 = 5; end31 = 10;
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
    loglog(dts,U_Err(i,:),'linestyle',linS{1},'color',C{i},'marker',Mar{i},'MarkerSize',ms)
    hold on
end
hold on
loglog(dts(st31:end31),Cof3(1)*dts(st31:end31).^Sl3(1),'--','color',Cref{1})
legend('','Location', 'best')
switch eqn
    case 'LinAdv'
        legend('(5,3,3)','(5,3,3)','(12,3,3)',sprintf('Slope %d',Sl3(1)),'NumColumns',1)
    case 'Burgers'
        legend('(5,3,3)','(5,3,3)','(12,3,3)',sprintf('Slope %d',Sl3(1)),'NumColumns',1)
end
xlim([dts(end),dts(1)])
ylim([lim3(1),lim3(2)])
xlabel('\Delta t');
ylabel('Error');
set(gca,'FontSize',fs)

% Save in eps
% figure_name = sprintf('Figures/%s_ConvPlotOrder_%d',eqn,P(1));
% print(gcf,figure_name,'-depsc','-r100','-vector')
%-------------------------------%
