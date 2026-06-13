function dydt=Raman_transition_new(t,y,gravity,tstart,H_B,pre,params,w_v)
% y1= [y(1:8)';horzcat(zeros(1,1),y(9:15)');horzcat(zeros(1,2),y(16:21)');
%     horzcat(zeros(1,3),y(22:26)');horzcat(zeros(1,4),y(27:30)');horzcat(zeros(1,5),y(31:33)');
%     horzcat(zeros(1,6),y(34:35)');horzcat(zeros(1,7),y(36)')];
% y2=conj(y1');
% y2=y2-diag(diag(y2));

% y=y1+y2;
v=y(4:6);
r=y(1:3);
rho=reshape(y(7:end),8,8);
tau1 =12e-6;
t=t-tstart;
A1 = 0.5*(1 + erf(t/(sqrt(2)*tau1)));
Aoff =1;

%
k1= params.beam1.khat*2*pi/780.24e-9;
k2= params.beam2.khat*2*pi/780.24e-9;
% 把F=1视为基线
H_detuning_unit = diag([1,1,1, 1,1,0,0,0]);

% --- 高斯强度与梯度（两束光） ---
[G1, dlnG1] = gaussian_profile_xy(r, params.beam1);
[G2, dlnG2] = gaussian_profile_xy(r, params.beam2);

% 局域拉比（Ω∝sqrt(I)）
% --- 局域拉比（两种方式二选一） ---
% 开关：params.use_I2 = 0(默认，直接用Omega2R) / 1(用光强重算)
%      （如需同样控制 Omega1，可设置 params.use_I1）
% 先算高斯包络对应的局域强度（W/m^2）
I1_loc = ( params.Itotal)/(1+params.I_ratio);   % 若未设置 I1_0 则为 0
I2_loc = (params.Itotal* params.I_ratio)/(1+params.I_ratio) ;
c0   = 2.99792458e8;
eps0 = 8.8541878128e-12;
hbar = 1.054571817e-34;
d    =  3.58e-29;
E1 = sqrt(2*I1_loc/(c0*eps0))*G1*A1*Aoff ;
E2 = sqrt(2*I2_loc/(c0*eps0))*G2*A1*Aoff ;

% Om1（可选，和 Om2 同逻辑；不想用就保留你原来的写法）
if isfield(params,'use_I1') && params.use_I1
    
    Om1 = d * E1 / hbar;%（单位 rad/s per sqrt(W/m^2)）');*A1
    Om2 =d * E2/ hbar;
else
    Om1 = params.Omega1R * sqrt(G1);  % 你原来的：按包络缩放
    Om2 = params.Omega1R * sqrt(G2);
end

% --- 用预计算基矩阵拼装局域 H_AC ---
H_AC = (abs(Om1)^2) * pre.M1 + (abs(Om2)^2) * pre.M2 ...
     + (conj(Om1)*Om2)*pre.M12* exp(-1i * params.phi_eff) + (Om1*conj(Om2))*pre.M12'* exp(1i *params.phi_eff);

doppler_shift =dot(k1-k2, v)+dot((k1+k2*cos(0.1)),v);
% --- 【核心修改】计算包含补偿的总双光子失谐 ---
delta_total_for_atom = (w_v*2*pi - 0*abs(doppler_shift));%+((abs(w_v)+2e4)*2*pi - abs(doppler_shift))+((abs(w_v)-2e4)*2*pi - abs(doppler_shift))
%delta_total_for_atom = 6.770457390295067e+06;
H_detuning_term = delta_total_for_atom * H_detuning_unit;
H=H_B+ H_AC+H_detuning_term;
% 交叉项梯度： (conjΩ1Ω2)* 0.5*(∇ln G1+∇ln G2) M12 + h.c.
hbar = 1.054571817e-34;
dHdx = (abs(Om1)^2)*dlnG1(1)*pre.M1 + (abs(Om2)^2)*dlnG2(1)*pre.M2 ...
     + 0.5*(conj(Om1)*Om2)*(dlnG1(1)+dlnG2(1))*pre.M12 ...
     + 0.5*(Om1*conj(Om2))*(dlnG1(1)+dlnG2(1))*pre.M12';
dHdy = (abs(Om1)^2)*dlnG1(2)*pre.M1 + (abs(Om2)^2)*dlnG2(2)*pre.M2 ...
     + 0.5*(conj(Om1)*Om2)*(dlnG1(2)+dlnG2(2))*pre.M12 ...
     + 0.5*(Om1*conj(Om2))*(dlnG1(2)+dlnG2(2))*pre.M12';
dHdz = (abs(Om1)^2)*dlnG1(3)*pre.M1 + (abs(Om2)^2)*dlnG2(3)*pre.M2 ...
     + 0.5*(conj(Om1)*Om2)*(dlnG1(3)+dlnG2(3))*pre.M12 ...
     + 0.5*(Om1*conj(Om2))*(dlnG1(3)+dlnG2(3))*pre.M12';

Fx = -hbar * real( trace( rho * dHdx * 2 ) );
Fy = -hbar * real( trace( rho * dHdy * 2 ) );
Fz = -hbar * real( trace( rho * dHdz* 2 ) );
F_opt = [Fx;Fy;Fz];
% --- 速度方程：g + 光势力/质量 ---
mRb = 1.44316060e-25;                % 87Rb 质量 [kg]
vdot = gravity + F_opt/mRb;


dydt=[v; vdot ; reshape((H*rho-rho*H)/1i,64,1)];
end
function [G, gradLnG, mask_aperture] = gaussian_profile_xy(r, beam)
% GAUSSIAN_PROFILE  任意方向高斯光束的相对强度与 ∇ln G，并支持光阑截断
%
% 输入：
%   r    : [3xN] 或 [Nx3] 位置 (m)
%   beam : 结构体
%          .center [3x1]     光束中心，即光轴经过的一点
%          .w0     scalar    光腰半径 (m)，按 1/e^2 强度半径定义
%          .khat   [3x1]     单位波矢方向
%          .zR     可选      瑞利长度，提供则考虑 w(z) 变化
%          .aperture_radius   可选，光阑半径 (m)
%          .aperture_diameter 可选，光阑直径 (m)
%
% 输出：
%   G              : 1xN，相对强度 I/I0
%   gradLnG        : [3xN]，∇ln G
%   mask_aperture  : 1xN，true 表示在光阑内，false 表示被光阑挡住

    % ---- 形状处理 ----
    if size(r,1) ~= 3 && size(r,2) == 3
        r = r.';
    end
    assert(size(r,1) == 3, 'r 必须是 3xN 或 Nx3');

    % ---- 基本量 ----
    r0   = beam.center(:);
    khat = beam.khat(:);
    khat = khat / norm(khat);              % 防止传入的 khat 未归一化

    d    = r - r0;                         % 3xN，位移向量
    P    = eye(3) - khat*khat.';           % 投影到横向平面的矩阵

    r_perp = P*d;                          % 横向分量
    s      = khat.'*d;                     % 轴向坐标，1xN
    r2     = sum(r_perp.^2, 1);            % 横向径向距离平方

    w0 = beam.w0;

    % ---- 是否考虑轴向光斑变化 ----
    use_axial = isfield(beam,'zR') && ~isempty(beam.zR) && isfinite(beam.zR);

    if use_axial
        zR = beam.zR;
        w  = w0 * sqrt(1 + (s./zR).^2);    % w(s)
        amp = (w0./w).^2;                  % 轴向强度因子
    else
        w = w0 + 0*s;
        amp = ones(1, size(r,2));
    end

    % ---- 高斯强度 ----
    G = amp .* exp(-2*r2./(w.^2));

    % ---- ∇ln G ----
    if use_axial
        % ln G = 2 ln(w0/w) - 2 r_perp^2 / w^2
        W = w.^2;

        % 横向梯度
        term_perp = (-4./W) .* r_perp;

        % 轴向修正
        dlnw_ds = s ./ (zR^2 + s.^2);
        Q = r2;
        axial = -2*dlnw_ds + 4*Q.*(w0^2/zR^2).*s./(W.^2);

        gradLnG = term_perp + khat .* axial;
    else
        gradLnG = (-4/(w0^2)) * r_perp;
    end

    % ============================================================
    % 光阑截断：横向半径超过光阑半径的点，光强置零
    % ============================================================
    mask_aperture = true(1, size(r,2));

    if isfield(beam, 'aperture_radius') && ~isempty(beam.aperture_radius)
        R_ap = beam.aperture_radius;
        mask_aperture = (r2 <= R_ap^2);

    elseif isfield(beam, 'aperture_diameter') && ~isempty(beam.aperture_diameter)
        R_ap = beam.aperture_diameter / 2;
        mask_aperture = (r2 <= R_ap^2);
    end

    % 光阑外光强置零
    G(~mask_aperture) = 0;

    % 光阑外 ∇lnG 理论上无定义，这里为了后续力计算稳定，置零
    gradLnG(:, ~mask_aperture) = 0;
end