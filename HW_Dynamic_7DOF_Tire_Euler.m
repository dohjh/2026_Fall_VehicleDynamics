clc; clear;
close all;

% Vehicle design parameter
m=1321; Iz=2120; 
a=1.07; b=1.53; L=a+b;
w=1.48;   % wheel tread
Jw=9.93; % tire inertia 
Td=0;   % driving torque/ Brake torque 
R=0.3;   % tire radius
hg=0.7;  % cg height
g=9.81;



steer_max = 5*pi/180; % [rad]
Vx = 10*1000/3600;  % m/s

% coefficient for Magic formula -----------------------------------
% 1. Traction model
b0 = 1.65; b1 = -7.6118; b2 = 1122.6; b3 = -7.36e-3; b4 = 144.82; 
b5 = -7.6614e-2; b6 = -3.86e-3; b7 = 8.5055e-2; b8 = 7.5719e-2;
% 2. Cornering model
a0 = 1.6929; a1 = -55.2084; a2 = 1271.28; a3 = 1601.8; a4 = 6.4946; 
a5 = 4.7966e-3; a6 = -0.3875; a7 = 1.0; a8 = -4.5399e-2;
%------------------------------------------------------------------


% initial value -----------------------------------------------------------
ws1 = Vx/R; ws2=ws1; ws3=ws1; ws4=ws1; % initial wheel speed [rad/s] 

t_end = 10
h=0.005; t=(0:h:t_end); % time step

y=zeros(length(t),7);
y(1,:)=[Vx,0,0,ws1,ws2,ws3,ws4];  % [Vx,Vy,yaw rate, wheel speed 1,2,3,4] 

Fz=zeros(1,4); 
Fz(1:2)=(b*m*g)/(2*L);  % Front weight
Fz(3:4)=(a*m*g)/(2*L);  % Rear weight

Fzc(1,:)=[(b*m*g)/(2*L) (b*m*g)/(2*L) (a*m*g)/(2*L) (a*m*g)/(2*L)];
Fxc(1,:)=[0 0 0 0]; Fyc(1,:)=[0 0 0 0]; 
alphac(1,:)=[0 0 0 0];
%--------------------------------------------------------------------------

% Select Scenario (1: Single Lane Change, 2: Double Lane Change)
scenario = 1;

for i=1:length(t)-1
     
    Fzc(i,:) = Fz;  % 수직하중 기록 [N]
    Fz=0.001*Fz; %Fz[kN] applied to Magic formula
    
    % Tire Magic formula --------------------------------------------------
    Cx = b0;
    Dx = b1*Fz.^2+b2*Fz;
    Bx = (b3*Fz.^2+b4*Fz)./(exp(b5*Fz)*Cx.*Dx);
    Ex = b6*Fz.^2+b7*Fz+b8;
    
    Cy = a0;
    Dy = a1*Fz.^2+a2*Fz;
    By = a3*sin(2*atan(Fz/a4))./(Cy*Dy);
    Ey = a6*Fz+a7;
    %----------------------------------------------------------------------
    
   
    % 차량동역학 7DOF 모델 코딩하시오.
    
    % 현재 State Variable
    vx = y(i,1); vy = y(i,2); yaw = y(i,3);
    ws = y(i,4:7);

    if scenario == 1
        % 1 : Single Lane Change
        if t(i) >= 2 && t(i) < 4
            delta = (t(i) - 2) / 2 * steer_max;
        elseif t(i) >= 4
            delta = steer_max;
        else
            delta = 0;
        end
    else
        % 2 : Double Lane Change
        if t(i) >= 2 && t(i) < 4
            delta = (t(i) - 2) / 2 * steer_max;
        elseif t(i) >= 4 && t(i) < 5
            delta = steer_max;
        elseif t(i) >= 5 && t(i) < 7
            delta = steer_max - (t(i) - 5) / 2 * steer_max;
        else
            delta = 0;
        end
    end
    steer_angle(i) = delta;

    % 1. 타이어 중심속도 (V_ix, V_iy) 및 슬립각(alpha) 계산부
    v1x = vx - (w/2)*yaw; v1y = vy + a*yaw;
    v2x = vx + (w/2)*yaw; v2y = vy + a*yaw;
    v3x = vx - (w/2)*yaw; v3y = vy - b*yaw;
    v4x = vx + (w/2)*yaw; v4y = vy - b*yaw;
    
    alpha = zeros(1,4);
    alpha(1) = delta - atan2(v1y, max(v1x, 0.01));
    alpha(2) = delta - atan2(v2y, max(v2x, 0.01));
    alpha(3) = 0 - atan2(v3y, max(v3x, 0.01));
    alpha(4) = 0 - atan2(v4y, max(v4x, 0.01));
    alphac(i,:) = alpha;
    
    % 타이어 중심 실제 속도 크기 (U_i) 계산
    U1 = sqrt(v1x^2 + v1y^2);
    U2 = sqrt(v2x^2 + v2y^2);
    U3 = sqrt(v3x^2 + v3y^2);
    U4 = sqrt(v4x^2 + v4y^2);
    U_arr = [U1, U2, U3, U4];
    
    % Longitudinal Slip Ratio 계산 
    s = zeros(1,4);
    for j = 1:4
        vw = R * ws(j); 
        v_tire_x = U_arr(j) * cos(alpha(j)); 
        
        if vw >= v_tire_x 
            s(j) = (vw - v_tire_x) / max(abs(vw), 0.01);
        else 
            s(j) = (vw - v_tire_x) / max(abs(v_tire_x), 0.01);
        end
    end
    sigmac(i,:) = s;  % 종방향 슬립률 기록
    
    alpha_deg = alpha * 180 / pi;
    s_pct = s * 100;
    
    % 2. Longitudinal, Lateral Tire Force 계산
    Fx = Dx .* sin(Cx .* atan(Bx .* s_pct - Ex .* (Bx .* s_pct - atan(Bx .* s_pct))));
    Fy = Dy .* sin(Cy .* atan(By .* alpha_deg - Ey .* (By .* alpha_deg - atan(By .* alpha_deg))));
    Fxc(i,:) = Fx; Fyc(i,:) = Fy;

    % 3. 7DoF 운동방정식
    sum_FX = -(Fy(1) + Fy(2))*sin(delta) + (Fx(1) + Fx(2))*cos(delta) + Fx(3) + Fx(4);
    sum_FY = (Fy(1) + Fy(2))*cos(delta) + (Fx(1) + Fx(2))*sin(delta) + Fy(3) + Fy(4);
    sum_MZ = a*((Fy(1) + Fy(2))*cos(delta) + (Fx(1) + Fx(2))*sin(delta)) ...
             - b*(Fy(3) + Fy(4)) ...
             + (w/2)*((Fy(1) - Fy(2))*sin(delta) - (Fx(1) - Fx(2))*cos(delta) - (Fx(3) - Fx(4)));

    % 가속도, 요각가속도 계산
    dot_vx = sum_FX/m + yaw*vy;
    dot_vy = sum_FY/m - yaw*vx;
    dot_yaw = sum_MZ/Iz;

    % 가속도 기록(LT 계산하려고)
    ax_c(i) = dot_vx;
    ay_c(i) = dot_vy;

    % 휠 각가속도 계산
    dot_ws = zeros(1,4);
    for j = 1:4
        dot_ws(j) = (Td/4 - Fx(j)*R) / Jw;
    end

    % 4. 다음 스텝 계산을 위한 State 업데이트
    y(i+1, 1) = vx + h * dot_vx;
    y(i+1, 2) = vy + h * dot_vy;
    y(i+1, 3) = yaw + h * dot_yaw;
    y(i+1, 4:7) = ws + h * dot_ws;

    % 5. LT 적용해서 다음 스텝 F_z 업데이트
    mf = m * (b/L);
    mr = m * (a/L);
    
    N1 = (mf/2)*g - mf*ay_c(i)*(hg/w) - (m/2)*ax_c(i)*(hg/L);
    N2 = (mf/2)*g + mf*ay_c(i)*(hg/w) - (m/2)*ax_c(i)*(hg/L);
    N3 = (mr/2)*g - mr*ay_c(i)*(hg/w) + (m/2)*ax_c(i)*(hg/L);
    N4 = (mr/2)*g + mr*ay_c(i)*(hg/w) + (m/2)*ax_c(i)*(hg/L);
    
    Fz = max([N1, N2, N3, N4], 100);
end

% 끝 배열 인덱스 보정
steer_angle(length(t)) = steer_angle(end);
ax_c(length(t)) = ax_c(end);
ay_c(length(t)) = ay_c(end);
Fxc(length(t),:) = Fxc(end,:);
Fyc(length(t),:) = Fyc(end,:);
Fzc(length(t),:) = Fzc(end,:);
alphac(length(t),:) = alphac(end,:);
sigmac(length(t),:) = sigmac(end,:);

% 차량동역학 운동 plot 하시오.
scen_name = {'Single Lane Change', 'Double Lane Change'};
figure('Name', sprintf('7DOF Simulation - %s', scen_name{scenario}), 'Position', [100 100 800 600]);

subplot(3,2,1);
plot(t, steer_angle * 180/pi, 'LineWidth', 1.5);
title('Steer Angle [deg]'); xlabel('t [s]'); ylabel('deg'); grid on; ylim([0 6]);

subplot(3,2,2);
plot(t, y(:,3), 'LineWidth', 1.5);
title('Yaw Rate [rad/s]'); xlabel('t [s]'); ylabel('rad/s'); grid on;

subplot(3,2,3);
plot(t, y(:,1) * 3600/1000, 'LineWidth', 1.5);
title('Vx [km/h]'); xlabel('t [s]'); ylabel('km/h'); grid on;

subplot(3,2,4);
plot(t, y(:,2) * 3600/1000, 'LineWidth', 1.5);
title('Vy [km/h]'); xlabel('t [s]'); ylabel('km/h'); grid on;

subplot(3,2,5);
plot(t, ax_c, 'LineWidth', 1.5);
title('ax [m/s^2]'); xlabel('t [s]'); ylabel('a_x'); grid on;

subplot(3,2,6);
plot(t, ay_c, 'LineWidth', 1.5);
title('ay [m/s^2]'); xlabel('t [s]'); ylabel('a_y'); grid on;

% 타이어 상태 plot 추가-> 보고서용(슬립각, 횡력, 종력, 수직하중)
tire_lbl = {'1 (FL)', '2 (FR)', '3 (RL)', '4 (RR)'};
figure('Name', sprintf('Tire States - %s', scen_name{scenario}), 'Position', [950 100 800 550]);

subplot(2,2,1);
plot(t, alphac * 180/pi, 'LineWidth', 1.3);
title('Slip Angle \alpha_i [deg]'); xlabel('t [s]'); ylabel('deg'); grid on;
legend(tire_lbl, 'Location', 'best');

subplot(2,2,2);
plot(t, Fyc, 'LineWidth', 1.3);
title('Lateral Force F_i [N]'); xlabel('t [s]'); ylabel('N'); grid on;

subplot(2,2,3);
plot(t, Fxc, 'LineWidth', 1.3);
title('Longitudinal Force D_i [N]'); xlabel('t [s]'); ylabel('N'); grid on;

subplot(2,2,4);
plot(t, Fzc, 'LineWidth', 1.3);
title('Normal Load N_i [N]'); xlabel('t [s]'); ylabel('N'); grid on;
