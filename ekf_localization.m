function [] = ekf_localization()
 
% Homework for ekf localization
% Modified by YH on 09/09/2019, thanks to the original open source
% Any questions please contact: zjuyinhuan@gmail.com

    close all;
    clear all;

    disp('EKF Start!')

    time = 0;
    global endTime; % [sec]
    endTime = 60;
    global dt;
    dt = 0.1; % [sec]

    removeStep = 5;

    nSteps = ceil((endTime - time)/dt);

    estimation.time=[];
    estimation.u=[];
    estimation.GPS=[];
    estimation.xOdom=[];
    estimation.xEkf=[];
    estimation.xTruth=[];

    % State Vector [x y yaw]'
    xEkf=[0 0 0]';
    PxEkf = eye(3);

    % Ground True State
    xTruth=xEkf;

    % Odometry Only
    xOdom=xTruth;

    % Observation vector [x y yaw]'
    z=[0 0 0]';

    % Simulation parameter
    global noiseQ
    noiseQ = diag([0.1 0 degreeToRadian(10)]).^2; %[Vx Vy yawrate]

    global noiseR
    noiseR = diag([0.5 0.5 degreeToRadian(5)]).^2;%[x y yaw]
    
    % Covariance Matrix for motion
    convQ=eye(3);

    % Covariance Matrix for observation
    convR=noiseR;

    % Other Intial
    % ?

    % Main loop
    for i=1 : nSteps
        time = time + dt;
        % Input
        u=robotControl(time);
        % Observation
        [z,xTruth,xOdom,u]=prepare(xTruth, xOdom, u);

        % ------ Kalman Filter --------
        % Predict
        % ?
        x_p = doMotion(xEkf,u);
        G_t = jacobF(xEkf,u);
        convQ = G_t * convR *G_t' + noiseQ;
        % Update
        % xEkf=?
        H = jacobH(x_p);
        S = H*convQ*H'+noiseR;
        K = convQ*H'*inv(S);
        xEkf = x_p +K*(z - doObservation(z,x_p));
        convR = (eye(3)-K*H)*convQ;


        % -----------------------------

        % Simulation estimation
        estimation.time=[estimation.time; time];
        estimation.xTruth=[estimation.xTruth; xTruth'];
        estimation.xOdom=[estimation.xOdom; xOdom'];
        estimation.xEkf=[estimation.xEkf;xEkf'];
        estimation.GPS=[estimation.GPS; z'];
        estimation.u=[estimation.u; u'];

        % Plot in real time
        % Animation (remove some flames)
        if rem(i,removeStep)==0
            %hold off;
            plot(estimation.GPS(:,1),estimation.GPS(:,2),'*m', 'MarkerSize', 5);hold on;
            plot(estimation.xOdom(:,1),estimation.xOdom(:,2),'.k', 'MarkerSize', 10);hold on;
            plot(estimation.xEkf(:,1),estimation.xEkf(:,2),'.r','MarkerSize', 10);hold on;
            plot(estimation.xTruth(:,1),estimation.xTruth(:,2),'.b', 'MarkerSize', 10);hold on;
            axis equal;
            grid on;
            drawnow;
            %movcount=movcount+1;
            %mov(movcount) = getframe(gcf);
        end 
    end
    close
    
    finalPlot(estimation);
 
end

% control
function u = robotControl(time)
    global endTime;

    T = 10; % sec
    Vx = 1.0; % m/s
    Vy = 0.2; % m/s
    yawrate = 5; % deg/s
    
    % half
    if time > (endTime/2)
        yawrate = -5;
    end
    
    u =[ Vx*(1-exp(-time/T)) Vy*(1-exp(-time/T)) degreeToRadian(yawrate)*(1-exp(-time/T))]';
    
end

% all observations for 
function [z, xTruth, xOdom, u] = prepare(xTruth, xOdom, u)
    global noiseQ;
    global noiseR;

    % Ground Truth
    xTruth=doMotion(xTruth, u);
    % add Motion Noises
    u=u+noiseQ*randn(3,1);
    % Odometry Only
    xOdom=doMotion(xOdom, u);
    % add Observation Noises
    z=xTruth+noiseR*randn(3,1);
end


% Motion Model
function x = doMotion(x, u)
    global dt;
    %?
    vel = sqrt(u(1)*u(1)+u(2)*u(2));
    vel_w = vel / u(3);
    temp_1 = -vel_w*sin(x(3))+vel_w*sin(x(3)+u(3)*dt);
    temp_2 = vel_w*cos(x(3))-vel_w*cos(x(3)+u(3)*dt);
    temp_3 = dt*u(3);
    temp = [temp_1 temp_2 temp_3]';
    x = x + temp;
end

% Jacobian of Motion Model
function jF = jacobF(x, u)
    global dt;
    vel = sqrt(u(1)*u(1)+u(2)*u(2));
    vel_w2 = vel / u(3);
    %dist = sqrt(u(1)*u(1) + u(2)*u(2))*dt;
    jF = [1 0 -vel_w2*cos(x(3))+vel_w2*cos(x(3)+u(3)*dt)
          0 1 -vel_w2*sin(x(3))+vel_w2*sin(x(3)+u(3)*dt)
          0 0 1  ];
    %?
end

%Observation Model
function x = doObservation(z, xPred)
    m = [1 0 0 
        0 1 0
        0  0 1];
    x=m*xPred;
    %?
 end

%Jacobian of Observation Model
function jH = jacobH(x)
    %?
    jH = [1 0 0
         0 1 0 
         0 0 1 ];
end

% finally plot the results
function []=finalPlot(estimation)
    figure;
    
    plot(estimation.GPS(:,1),estimation.GPS(:,2),'*m', 'MarkerSize', 5);hold on;
    plot(estimation.xOdom(:,1), estimation.xOdom(:,2),'.k','MarkerSize', 10); hold on;
    plot(estimation.xEkf(:,1), estimation.xEkf(:,2),'.r','MarkerSize', 10); hold on;
    plot(estimation.xTruth(:,1), estimation.xTruth(:,2),'.b','MarkerSize', 10); hold on;
    legend('GPS Observations','Odometry Only','EKF Localization', 'Ground Truth');

    xlabel('X (meter)', 'fontsize', 12);
    ylabel('Y (meter)', 'fontsize', 12);
    grid on;
    axis equal;
    
    % calculate error
    % ?
    num = size(estimation.xTruth,1);
    fprintf("终点站的真值：x = %f,y = %f\n",estimation.xTruth(num,1),estimation.xTruth(num,2));
    fprintf("终点站的EKF：x = %f,y = %f\n",estimation.xEkf(num,1),estimation.xEkf(num,2));
    fprintf("终点站的里程计Odom：x = %f,y = %f\n",estimation.xOdom(num,1),estimation.xOdom(num,2));
    miao = 1;error_Odom = 0;error_Ekf = 0;
    
    while miao <= num
        error_Odom = error_Odom + sqrt((estimation.xOdom(miao,1)-estimation.xTruth(miao,1)).^2 + (estimation.xOdom(miao,2)-estimation.xTruth(miao,2)).^2);
        error_Ekf = error_Ekf +sqrt((estimation.xEkf(miao,1)-estimation.xTruth(miao,1)).^2 + (estimation.xEkf(miao,2)-estimation.xTruth(miao,2)).^2);
        miao = miao + 1;
    end
    disp(["num=",num]);
    fprintf("Odom累计平均误差：%f\n",error_Odom/num);
    fprintf("Ekf累计平均误差：%f\n",error_Ekf/num);

        


end

function radian = degreeToRadian(degree)
    radian = degree/180*pi;
end