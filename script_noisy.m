%% Generates Figures 4 and 5
close all; clear; clc

% save_figs = true;
save_figs = false;

path_name = '../../Figures/'; % save figures here
Nsims = 10000; % number of realisations
colors = parula(4); % colors for plots

%% Material and laser flash parameters
Qinf = 7000; % total heat absorbed
l0 = 0.001; % inner surface at x = l0
l1 = 0.003; % outer surface at x = l1
k = 222; % thermal conductivity
rho = 2700; % density
cp = 896; % specific heat capacity
N = 1000; % number of temperature rise values (excluding t = 0)
tN = 0.1; % end time
beta = 0.001; % exponential pulse parameter (peak occurs at t = beta)
q = @(t) Qinf*t.*exp(-t/beta)/beta^2; % exponential pulse
dt = tN/N; % time step duration
t = (0:dt:tN)'; % discrete times
alpha = k/(rho*cp); % target value of thermal diffusivity
sigma = [0.02,0.05]; % standard deviation values for noise

%% Finite volume method parameters
Nx = 501; % number of nodes
h = (l1-l0)/(Nx-1); % node spacing
x = linspace(l0,l1,Nx); % location of nodes
xw(1) = x(1); xw(2:Nx) = (x(1:Nx-1)+x(2:Nx))/2; % west boundaries
xe(1:Nx-1) = xw(2:Nx); xe(Nx) = x(Nx); % east boundaries
T0 = zeros(1,Nx); % initial temperature rise
AbsTol = 1e-12; % absolute error tolerence
RelTol = 1e-12; % relative error tolerence

% Sparsity pattern of Jacobian
e = ones(Nx,1); JPat = spdiags([e e e],-1:1,Nx,Nx);
options = odeset('RelTol',RelTol,'AbsTol',AbsTol,'JPattern',JPat);

%% Generates Figures 4(a),4(b),4(c),4(d)
alpha_tilde = zeros(Nsims,length(sigma),3);
epsilon = zeros(2,Nsims,length(sigma),3);
configurations = {'outward','inward'};
rng default; rng(1); rand_matrix = randn(N+1,Nsims);

for kk = 1:2
    
    configuration = configurations{kk};
    
    for j = 1:length(sigma)
        figure;
        scrz = get(gcf,'OuterPosition');
        set(gcf,'OuterPosition',[scrz(1:2) scrz(3)*1.0 scrz(4)],'Color','w')
        
        noise = sigma(j)*rand_matrix;
        
        for d = 1:3
            
            % Finite pulse (FVM)
            [~,T] = ode15s(@(t,T) Gfunc(t,T,d,alpha,l0,l1,h,x,xw,xe,rho,cp,Nx,configuration,q),t,T0,options);
            
            if isequal(configuration,'outward')
                Tdata = T(:,end);
                Tinf = d*l0^(d-1)*Qinf/(rho*cp*(l1^d-l0^d));
            elseif isequal(configuration,'inward')
                Tdata = T(:,1);
                Tinf = d*l1^(d-1)*Qinf/(rho*cp*(l1^d-l0^d));
            end
            
            for i = 1:Nsims
                Tdata_noisy = Tdata + noise(:,i);
                alpha_tilde(i,j,d) = thermal_diffusivity(d,t,Tdata_noisy,Tinf,l0,l1,beta);
                epsilon(kk,i,j,d) = (alpha - alpha_tilde(i,j,d))/alpha * 100; 
            end
            
            % experimental temperature rise history for last trial
            plot(t,Tdata_noisy,'-','LineWidth',1,'Color',colors(d,:));
            
            hold on
            % theorectical temperaure rise history for last estimated value of alpha
            [t,T] = ode15s(@(t,T) Gfunc(t,T,d,alpha_tilde(end,j,d),l0,l1,h,x,xw,xe,rho,cp,Nx,configuration,q),t,T0,options);
            if isequal(configuration,'outward')
                Tmodel = T(:,end);
            elseif isequal(configuration,'inward')
                Tmodel = T(:,1);
            end
            plot(t,Tmodel,'-','LineWidth',1,'Color','k')
            ex = -5; % exponent
            if isequal(configuration,'outward')
                text(0.2,0.9,['Outward [$\sigma = ',num2str(sigma(j),'%g'),'\,^{\circ}\mathrm{C}$]'],...
                'Color','k','Interpreter','LaTeX','FontSize',24,'Units','Normalized')
                yleg = (0.81-0.09*(d-1));
                plot([0.13,0.18]*tN,(yleg*3.7-0.3)*ones(1,2),'Color',colors(d,:),'LineWidth',4)
                text(0.2,yleg,['$d = ',num2str(d),'\!: \widetilde{\alpha} = ',num2str(alpha_tilde(end,j,d)/10^ex,'%1.4f'),...
                    ' \times 10^{',num2str(ex,'%i'),'}\,\mathrm{m}^{2}\mathrm{s}^{-1}$'],...
                    'Color','k','Interpreter','LaTeX','FontSize',24,'Units','Normalized')                
            elseif isequal(configuration,'inward')
                text(0.2,0.34,['Inward [$\sigma = ',num2str(sigma(j),'%g'),'\,^{\circ}\mathrm{C}$]'],...
                    'Color','k','Interpreter','LaTeX','FontSize',24,'Units','Normalized')
                yleg = (0.26-0.09*(d-1));
                plot([0.13,0.18]*tN,(yleg*3.7-0.3)*ones(1,2),'Color',colors(d,:),'LineWidth',4)
                text(0.2,yleg,['$d = ',num2str(d),'\!: \widetilde{\alpha} = ',num2str(alpha_tilde(end,j,d)/10^ex,'%1.4f'),...
                    ' \times 10^{',num2str(ex,'%i'),'}\,\mathrm{m}^{2}\mathrm{s}^{-1}$'],...
                    'Color','k','Interpreter','LaTeX','FontSize',24,'Units','Normalized')
            end
            ylims = get(gca,'YLim');
            xlim([0 tN]); ylim([-0.3 3.4]);
            set(gca,'FontSize',28,'TickLabelInterpreter','latex','XTick',[0,tN/2,tN],'Ytick',...
                [0,1.6,3.2])
            xl = xlabel('Time ($\mathrm{s}$)','Interpreter','LaTeX','FontSize',28);
            yl = ylabel('Temperature rise ($^{\circ}\mathrm{C}$)','Interpreter','LaTeX','FontSize',28);
            set(yl, 'Units','normalized');
            set(xl, 'Units','normalized');
            labels = {'(a)','(c)','(b)','(d)'};
            text(-0.185,-0.185,labels{j+2*isequal(configuration,'inward')},'FontSize',28, 'Units','normalized') 
            drawnow            
        end
        
        drawnow
        if save_figs
            pause(1);
            print(gcf,[path_name,'Tdata_','sigma',num2str(j),'_',configuration],'-depsc2')
        end
        
    end
    
end

%% Generates Figures 5(a),5(b),5(c),5(d)

for kk = 1:length(configurations)
    
    configuration = configurations{kk};
    
    for j = 1:length(sigma)
        
        figure;
        scrz = get(gcf,'OuterPosition');
        set(gcf,'OuterPosition',[scrz(1:2) scrz(3)*1.0 scrz(4)],'Color','w')
        set(gcf,'Renderer','Painters');
        
        % Smooth density profiles
        if isequal(configuration,'outward')
            [f,xi] = ksdensity(epsilon(kk,:,j,3)); plot(xi,f,'LineWidth',2), hold on
            [f,xi] = ksdensity(epsilon(kk,:,j,2)); plot(xi,f,'LineWidth',2)
            [f,xi] = ksdensity(epsilon(kk,:,j,1)); plot(xi,f,'LineWidth',2)
            [f,xi] = ksdensity(epsilon(kk,:,j,1)); area(xi,f,'LineWidth',2,'FaceColor',colors(1,:),'EdgeColor','none','FaceAlpha',0.25), hold on
            [f,xi] = ksdensity(epsilon(kk,:,j,2)); area(xi,f,'LineWidth',2,'FaceColor',colors(2,:),'EdgeColor','none','FaceAlpha',0.25)
            [f,xi] = ksdensity(epsilon(kk,:,j,3)); area(xi,f,'LineWidth',2,'FaceColor',colors(3,:),'EdgeColor','none','FaceAlpha',0.25)
            [f,xi] = ksdensity(epsilon(kk,:,j,1)); plot(xi,f,'LineWidth',2,'Color',colors(1,:)), hold on
            [f,xi] = ksdensity(epsilon(kk,:,j,2)); plot(xi,f,'LineWidth',2,'Color',colors(2,:))
            [f,xi] = ksdensity(epsilon(kk,:,j,3)); plot(xi,f,'LineWidth',2,'Color',colors(3,:))
            text(0.07,0.9,'Outward','Interpreter','LaTeX','FontSize',24,'Color','k','Units','normalized','HorizontalAlignment','left');
            text(0.07,0.81,['$\sigma = ',num2str(sigma(j),'%g'),'\,^{\circ}\mathrm{C}$'],'Interpreter','LaTeX','FontSize',24,...
                'Color','k','Units','normalized','HorizontalAlignment','left');
        elseif isequal(configuration,'inward')
            [f,xi] = ksdensity(epsilon(kk,:,j,1)); plot(xi,f,'LineWidth',2), hold on
            [f,xi] = ksdensity(epsilon(kk,:,j,2)); plot(xi,f,'LineWidth',2)
            [f,xi] = ksdensity(epsilon(kk,:,j,3)); plot(xi,f,'LineWidth',2)
            [f,xi] = ksdensity(epsilon(kk,:,j,3)); area(xi,f,'LineWidth',2,'FaceColor',colors(3,:),'EdgeColor','none','FaceAlpha',0.25), hold on
            [f,xi] = ksdensity(epsilon(kk,:,j,2)); area(xi,f,'LineWidth',2,'FaceColor',colors(2,:),'EdgeColor','none','FaceAlpha',0.25)
            [f,xi] = ksdensity(epsilon(kk,:,j,1)); area(xi,f,'LineWidth',2,'FaceColor',colors(1,:),'EdgeColor','none','FaceAlpha',0.25)
            [f,xi] = ksdensity(epsilon(kk,:,j,3)); plot(xi,f,'LineWidth',2,'Color',colors(3,:)), hold on
            [f,xi] = ksdensity(epsilon(kk,:,j,2)); plot(xi,f,'LineWidth',2,'Color',colors(2,:))
            [f,xi] = ksdensity(epsilon(kk,:,j,1)); plot(xi,f,'LineWidth',2,'Color',colors(1,:))
            text(0.07,0.9,'Inward','Interpreter','LaTeX','FontSize',24,'Color','k','Units','normalized','HorizontalAlignment','left');
            text(0.07,0.81,['$\sigma = ',num2str(sigma(j),'%g'),'\,^{\circ}\mathrm{C}$'],'Interpreter','LaTeX','FontSize',24,...
                'Color','k','Units','normalized','HorizontalAlignment','left');
        end
        
        % 0.5% and 99.5% quantiles
        q1 = round(quantile(epsilon(kk,:,j,1),[0.005,0.995]),2);
        q2 = round(quantile(epsilon(kk,:,j,2),[0.005,0.995]),2);
        q3 = round(quantile(epsilon(kk,:,j,3),[0.005,0.995]),2);
        
        yleg = 0.9;
        plot([0.6,0.65]*20-10,(yleg*1.3)*ones(1,2),'Color',colors(1,:),'LineWidth',4)
        text(0.67,0.9,'$d = 1$','Interpreter','LaTeX','FontSize',24,'Color','k','Units','normalized','HorizontalAlignment','left');
        text(0.6,0.81,['[',num2str(q1(1),'%.2f'),', ',num2str(q1(2),'%.2f'),']\%'],'Interpreter','LaTeX','FontSize',24,...
            'Color','k','Units','normalized','HorizontalAlignment','left');
        yleg = 0.69;
        plot([0.6,0.65]*20-10,(yleg*1.3)*ones(1,2),'Color',colors(2,:),'LineWidth',4)
        text(0.67,0.69,'$d = 2$','Interpreter','LaTeX','FontSize',24,'Color','k','Units','normalized','HorizontalAlignment','left');
        text(0.6,0.6,['[',num2str(q2(1),'%.2f'),', ',num2str(q2(2),'%.2f'),']\%'],'Interpreter','LaTeX','FontSize',24,...
            'Color','k','Units','normalized','HorizontalAlignment','left');
        yleg = 0.48;
        plot([0.6,0.65]*20-10,(yleg*1.3)*ones(1,2),'Color',colors(3,:),'LineWidth',4)
        text(0.67,0.48,'$d = 3$','Interpreter','LaTeX','FontSize',24,'Color','k','Units','normalized','HorizontalAlignment','left');
        text(0.6,0.39,['[',num2str(q3(1),'%.2f'),', ',num2str(q3(2),'%.2f'),']\%'],'Interpreter','LaTeX','FontSize',24,...
            'Color','k','Units','normalized','HorizontalAlignment','left');
        ylims = get(gca,'YLim');
        set(gca,'FontSize',28,'TickLabelInterpreter','latex','XTick',[-10,-5,0,5,10],'Ytick',[0,0.6,1.2])
        xlabel('$\varepsilon$ (\%)','Interpreter','LaTeX','FontSize',28)
        ylabel('Density','Interpreter','LaTeX','FontSize',28)
        labels = {'(a)','(c)','(b)','(d)'};
        text(-0.185,-0.185,labels{j+2*isequal(configuration,'inward')},'FontSize',28, 'Units','normalized')
        xlim([-10 10]); ylim([0 1.3]);
        drawnow
        
        if save_figs
            pause(1);
            print(gcf,[path_name,'Errors_','sigma',num2str(j),'_',configuration],'-depsc2')
        end
        
        pause(0.1)
        
    end
    
end

