% Example generation:
%     dt = .01; T = 0:dt:10;   % simulation time vector
    dt = .005; T = 0:dt:1;   % simulation time vector
%     alpha = 0.9;            % initial guess of alpha
%     beta = 0.005;           % initial guess of beta
    XR = zeros(length(T),2); XM = zeros(length(T),2);
    Trajectory_origin = [1 1];
%     rozvoj_amplitudy=1;rozvoj_faze=1;
    rozvoj_amplitudy=10;rozvoj_faze=1;
    err_mean_x = 0.0;
    err_mean_y = 0.0;
    std_x = 0.2;
    std_y = 0.2;
    for t = 1:length(T)
%         xm = [5*cos(2*pi*T(t)) 5*sin(2*pi*T(t))]';  % true position (A circle)
%         xm = [(1*T(t)+0.01*randn)*cos(2*pi*T(t)) (1*T(t)+0.05*randn)*sin(2*pi*T(t))]';  % true position
%         xm = [(Trajectory_origin(1)+1*T(t))*cos(2*pi*T(t)) (Trajectory_origin(2)+1*T(t))*sin(2*pi*T(t))]';  % true position (A spiral)
        xm = [(Trajectory_origin(1)+rozvoj_amplitudy*T(t))*cos(2*pi*T(t)/rozvoj_faze) (Trajectory_origin(2)+rozvoj_amplitudy*T(t))*sin(2*pi*T(t)/rozvoj_faze)]';  % true position (A spiral)
%         xm = [(Trajectory_origin(1)+rozvoj_amplitudy*T(t)) (Trajectory_origin(2)+rozvoj_amplitudy*T(t))*1.5]';  % true position (linear movement)
        XR(t,:) = xm;
        xm(1,1) = xm(1,1) + err_mean_x + std_x .* randn;       % x with error mean x and sd x
        xm(2,1) = xm(2,1) + err_mean_y + std_y .* randn;       % y with error mean y and sd y
        XM(t,:) = xm;
    end
%     save('Sem10Data_001.mat','dt','T','XR','XM','err_mean_x','err_mean_y','std_x','std_y');
%     save('Sem10Data_002.mat','dt','T','XR','XM','err_mean_x','err_mean_y','std_x','std_y');
%     save('Sem10Data_003.mat','dt','T','XR','XM','err_mean_x','err_mean_y','std_x','std_y');

%% Processing
% load('Sem10Data_001_linear.mat');
% load('Sem10Data_002_circular.mat');
load('Sem10Data_003_spiral2.mat');
% load('Sem10Data_003_spiral1.mat');
% load('Sem10Data_003.mat');
 XK = zeros(length(T),2); VK = zeros(length(T),2); RK = zeros(length(T),2);
%     xk = [0 0]';            % initial state (position)
    xk = XM(1,:)';            % initial state (position)
    vk = [0 0]';            % initial dState/dt (velocity)
% good for linear Data_001
    alpha = 0.12;            % initial guess of alpha
    beta = 0.003;           % initial guess of beta
% good for linear Data_001
%     alpha = 0.2;            % initial guess of alpha
%     beta = 0.002;           % initial guess of beta   
    
    for t = 1:length(T)
%         xm = [cos(2*pi*T(t)) sin(2*pi*T(t))]';  % true position (A circle)
%         xm = [(1*T(t)+0.01*randn)*cos(2*pi*T(t)) (1*T(t)+0.05*randn)*sin(2*pi*T(t))]';  % true position
%         xm = [(Trajectory_origin(1)+1*T(t))*cos(2*pi*T(t)) (Trajectory_origin(2)+1*T(t))*sin(2*pi*T(t))]';  % true position (A spiral)
%         xm = [(Trajectory_origin(1)+rozvoj_amplitudy*T(t))*cos(2*pi*T(t)/rozvoj_faze) (Trajectory_origin(2)+rozvoj_amplitudy*T(t))*sin(2*pi*T(t)/rozvoj_faze)]';  % true position (A spiral)
        xm(1,1) = XM(t,1);
        xm(2,1) = XM(t,2);
        [xkp,vkp,rk] = alphaBetaFilter(xm, dt, xk, vk, alpha, beta);
        xk = xkp;
        vk = vkp;
        XK(t,:) = xkp; VK(t,:) = vkp; RK(t,:) = rk;
    end

%%
%     figure('units','pixels','Position',[0 0 1024 768]);
    figure(1);
    subplot(2,1,1);
%     plot(T,XR(:,1),T,XM(:,1),T,XK(:,1),T,VK(:,1));
%     xlabel('Time (s)'); title('X Pos, Vel and Acc');
%     legend('Real','Measured','Estimated','Velocity');
    plot(T,XR(:,1),T,XM(:,1),T,XK(:,1));
    xlabel('Time (s)'); title('X Pos');
    legend('Real','Measured','Estimated');
    subplot(2,1,2);
%     plot(T,XR(:,2),T,XM(:,2),T,XK(:,2),T,VK(:,2));
%     xlabel('Time (s)'); title('Y Pos, Vel and Acc');
%     legend('Real','Measured','Estimated','Velocity');
    plot(T,XR(:,2),T,XM(:,2),T,XK(:,2));
    xlabel('Time (s)'); title('Y Pos');
    legend('Real','Measured','Estimated');

    figure(2)
    plot(XR(:,1),XR(:,2))
    hold on
    plot(XM(:,1),XM(:,2),'+')
    plot(XK(:,1),XK(:,2),'-*')
    hold off
    
    chyba=(XR-XK);
    chyba2=chyba.*chyba;
    chyba3=chyba2*[1;1];
    chyba4=sqrt(chyba3);
    mchyba=(XR-XM);
    mchyba2=mchyba.*mchyba;
    mchyba3=mchyba2*[1;1];
    mchyba4=sqrt(mchyba3);
    figure(3)
    plot(T,mchyba4,'r')
    hold on
    plot(T,chyba4,'y')
    hold off
    legend('RMS meas-real','RMS est-real')