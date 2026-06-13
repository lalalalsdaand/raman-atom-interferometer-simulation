function H = get_zeeman_hamiltonian(params, g, Rb_Quantum, B_vec)
% H  = get_zeeman_hamiltonian(params, g, Rb_Quantum, B_vec)
% 87Rb 5S1/2 (J=1/2, I=3/2) 在 |F,mF> 基（默认以 z 为量子化轴）的
% 超精细 + 线性塞曼 + 二阶塞曼（Löwdin 消元）哈密顿量（单位 rad/s）
%
% 额外支持：params.quant_axis = [qx,qy,qz] （量子化轴，默认 [0,0,1]）
% 若量子化轴不是 z，则先把实验室的 B_vec 旋转到该轴的坐标系中，再与 S/I 算符点乘。

    % ===== 常数 =====
    C.HBAR = 1.054571817e-34;     % J*s
    C.MU_B = 9.2740100783e-24;    % J/T

    % ===== 默认参数 =====
    if ~isfield(params,'B_unit') || isempty(params.B_unit), params.B_unit = 'G'; end
    if ~isfield(params,'delta_hfs_GHz') || isempty(params.delta_hfs_GHz)
        params.delta_hfs_GHz = 6.834682610;  % 87Rb ground-state HFS (GHz)
    end
    if ~isfield(params,'quant_axis') || isempty(params.quant_axis)
        params.quant_axis = [0;0;1];
    end
    if ~isfield(Rb_Quantum,'gJ') || isempty(Rb_Quantum.gJ), Rb_Quantum.gJ = 2.00233113; end
    if ~isfield(Rb_Quantum,'gI') || isempty(Rb_Quantum.gI), Rb_Quantum.gI = -0.0009951414; end
    if ~isfield(Rb_Quantum,'I')  || isempty(Rb_Quantum.I),  Rb_Quantum.I  = 3/2; end
    if ~isfield(Rb_Quantum,'J')  || isempty(Rb_Quantum.J),  Rb_Quantum.J  = 1/2; end

    % ===== 单位换算（Gauss/Tesla） =====
    B_vec = B_vec(:);
    switch upper(params.B_unit)
        case 'G',  B_T_lab = 1e-4 * B_vec;   % Gauss -> Tesla
        case 'T',  B_T_lab = B_vec;
        otherwise, error('未知 B_unit=%s，仅支持 ''G'' 或 ''T''。', params.B_unit);
    end

    % ===== 把实验室 B 旋到“量子化轴坐标系” =====
    % 量子化轴 q̂ 作为新坐标系的 z'，构造将 lab 坐标 -> {x',y',z'=q̂} 的旋转 R
    qhat = params.quant_axis(:);  qhat = qhat / norm(qhat);
    R = lab_to_quant_rotation(qhat);     % 3x3
    B_T = R * B_T_lab;                   % 用这个 B_T = [Bx',By',Bz'] 去点乘 S/I

    % ===== 角频率常数 =====
    muB_over_hbar = C.MU_B / C.HBAR;     % (rad/s)/T
    gJ = Rb_Quantum.gJ;  gI = Rb_Quantum.gI;

    % ===== 超精细常数（rad/s）=====
    delta_hfs_Hz = params.delta_hfs_GHz * 1e9;     % Hz
    % F=2 与 F=1 的角频率差 Δω = 2π Δν
    Delta = 2*pi * delta_hfs_Hz;                   % rad/s
    % 用 A/4 * D 形式：Δω = 2A_w  ->  A_w = Δω/2
    A_w = Delta/2;                                  % rad/s
    % D 的对角元：F=2 -> +3，F=1 -> -5
    D_diag = zeros(g.N,1); D_diag(g.F==2)=+3; D_diag(g.F==1)=-5;
    H0 = (A_w/4) * diag(D_diag);

    % ===== 生成算符 =====
    F_ops = get_F_operators(g);    Fx = F_ops{1}; Fy = F_ops{2}; Fz = F_ops{3};
    S_ops = get_spin_operators(g, Rb_Quantum.I, Rb_Quantum.J);
    Sx = S_ops.S{1}; Sy = S_ops.S{2}; Sz = S_ops.S{3};
    Ix = Fx - Sx; Iy = Fy - Sy; Iz = Fz - Sz;

    % ===== Zeeman 相互作用 V = (μB/ħ) [ gJ S + gI I ] · B' =====
    SB = B_T(1)*Sx + B_T(2)*Sy + B_T(3)*Sz;
    IB = B_T(1)*Ix + B_T(2)*Iy + B_T(3)*Iz;
    V  = muB_over_hbar * ( gJ*SB + gI*IB );     % rad/s

    % ===== 投影器与块分解 =====
    P2 = diag(double(g.F==2));   % 5x5
    P1 = diag(double(g.F==1));   % 3x3
    V_parallel = P2*V*P2+ P1*V*P1;         % 块内（不混 F）
    V_perp     = V - V_parallel;            % 跨块（F 混合）

    % ===== 二阶有效哈密顿量（Löwdin/Schrieffer–Wolff）=====
    % H^(2) = (1/Δ) [ P2 V⊥ P1 V⊥ P2  -  P1 V⊥ P2 V⊥ P1 ]
    H2 = (1/Delta) * ( P2*V_perp*P1*V_perp*P2 - P1*V_perp*P2*V_perp*P1 );

    % ===== 总哈密顿量 =====
    H = V_parallel + H2;    % rad/s

    % 厄米性检查
    nh = norm(H - H','fro');
    if nh > 1e-9*max(1,norm(H,'fro'))
        warning('H 非厄米：||H-H^†||=%.3g', nh);
    end
end

function R = lab_to_quant_rotation(qhat)
% 返回把“实验室坐标”旋到“量子化轴坐标系{x',y',z'=q̂}”的 3x3 旋转矩阵
    z = [0;0;1];
    if norm(cross(z,qhat)) < 1e-12
        % 与 z 同向或反向
        if dot(z,qhat) > 0, R = eye(3); else, R = diag([1,-1,-1]); end
        return;
    end
    v = cross(z,qhat);  s = norm(v);  c = dot(z,qhat);
    vx = [   0, -v(3),  v(2);
           v(3),   0 , -v(1);
          -v(2), v(1),   0  ];
    R = eye(3) + vx + vx*vx * ((1-c)/(s^2));  % Rodrigues
end
