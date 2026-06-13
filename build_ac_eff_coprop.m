function [dHg_plus, dHg_minus, dHg_0,parts] = build_ac_eff_coprop(g, C, params)
% BUILD_AC_EFF 生成 Rb87 D2 的有效(AC + Raman)哈密顿量分块
%
% 输入:
%   g      : 结构体，至少可含 g.N, g.F
%   C      : 原子常数与能级参数
%   params : 结构体，需包含
%            .delta1, .delta2
%            .e1_vec, .e2_vec
%            .Omega1R, .Omega2R
%            .phi_eff
%            .quant_axis
%            .include_mutual_ac
%
% 输出:
%   dHg_minus : e1 在 F=1 主作用 AC
%   dHg_plus  : e2 在 F=2 主作用 AC
%   dHg_0     : 互扰 AC + Raman 交叉
%   parts     : 分项结构体

    if nargin < 3, params = struct; end
    if ~isfield(params,'quant_axis'),        params.quant_axis = [0;0;1]; end
    if ~isfield(params,'include_mutual_ac'), params.include_mutual_ac = true; end
    if ~isfield(params,'phi_eff'),           params.phi_eff = 0; end

    N  = getfield_or(g,'N',8);
    Fv = getfield_or(g,'F',[2 2 2 2 2 1 1 1].');
    Pa = diag(Fv==1);   % F=1 projector
    Pb = diag(Fv==2);   % F=2 projector

    e1_in = params.e1_vec;
    e2_in = params.e2_vec;
    qa    = params.quant_axis;
    Om1   = params.Omega1R;
    Om2   = params.Omega2R;
    phi   = params.phi_eff;
    w_hfs = C.W_HFS;

    % ---------- detuning structs ----------
    det12   = struct('delta1', params.delta1,           'delta2', params.delta2);
    det1_F1 = struct('delta1', params.delta1,           'delta2', params.delta1);
    det1_F2 = struct('delta1', params.delta1 +w_hfs,   'delta2', params.delta1 +w_hfs );
    det2_F2 = struct('delta1', params.delta2,           'delta2', params.delta2);
    det2_F1 = struct('delta1', params.delta2 - w_hfs,   'delta2', params.delta2 - w_hfs);

    % ---------- tensors ----------
%     opt12 = struct('mode','raman','rule','avg');
%     [T12, ~] = dipole_T_FF_auto(e1_in, e2_in, qa, C, det12, opt12);
% 
%     opt11 = struct('mode','raman','rule','delta1');
%     [T11_F1, ~] = dipole_T_FF_auto(e1_in, e1_in, qa, C, det1_F1, opt11);  %F1光F1跃迁
%     [T11_F2, ~] = dipole_T_FF_auto(e1_in, e1_in, qa, C, det1_F2, opt11);  %F2光F1跃迁
%     opt11 = struct('mode','raman','rule','delta1');
%     [T2_on_F1, ~] = dipole_T_FF_auto(e2_in, e2_in, qa, C, det1_F1, opt11);
%     [T2_on_F2, ~] = dipole_T_FF_auto(e2_in, e2_in, qa, C, det1_F2, opt11);
% 
%     opt22 = struct('mode','raman','rule','delta2');
%     [T22_F2, ~] = dipole_T_FF_auto(e2_in, e2_in, qa, C, det2_F2, opt22);
%     [T22_F1, ~] = dipole_T_FF_auto(e2_in, e2_in, qa, C, det2_F1, opt22);
%     [T1_on_F2, ~] = dipole_T_FF_auto(e1_in, e1_in, qa, C, det2_F2, opt22);
%     [T1_on_F1, ~] = dipole_T_FF_auto(e1_in, e1_in, qa, C, det2_F2, opt22);

    [D1, P1] = build_D_from_pol_auto(e1_in, qa);
    [D2, P2] = build_D_from_pol_auto(e2_in,  qa);
    Fp = P1.Fp;   % 16x1
   % ---------- Raman 非对角 ----------
    [W12, info12] = build_W_from_detune(Fp, C, det12, 'avg');
    T12 = build_T_from_D(D1, D2, W12);
    T21 = build_T_from_D(D2, D1, W12);

    % ---------- beam1 的单束光移 ----------
    [W1_F1, info1F1] = build_W_from_detune(Fp, C, det1_F1, 'delta1');
    [W1_F2, info1F2] = build_W_from_detune(Fp, C, det1_F2, 'delta1');

    T1_on_F1 = build_T_from_D(D1, D1, W1_F1);
    T1_on_F2 = build_T_from_D(D1, D1, W1_F2);

    [W2_F1, info2F1] = build_W_from_detune(Fp, C, det2_F1, 'delta2');
    [W2_F2, info2F2] = build_W_from_detune(Fp, C, det2_F2, 'delta2');

    T2_on_F1 = build_T_from_D(D2, D2, W2_F1);
    T2_on_F2 = build_T_from_D(D2, D2, W2_F2);


%     % ---------- beam2 的单束光移 ----------
%     [W2_F1, info2F1] = build_W_from_detune(Fp, C, det2_F1, 'delta2');
%     [W2_F2, info2F2] = build_W_from_detune(Fp, C, det2_F2, 'delta2');
% 
%     T2_on_F1 = build_T_from_D(D2, D2, W2_F1);
%     T2_on_F2 = build_T_from_D(D2, D2, W2_F2);
% 
%      [W1_F1_1, info1F1] = build_W_from_detune(Fp, C, det1_F1, 'delta1'); %对于F1光
%     [W1_F2_1, info1F2] = build_W_from_detune(Fp, C, det1_F2, 'delta1');
% 
% 
%     T1_on_F1_1 = build_T_from_D(D2, D2, W1_F1_1);
%     T1_on_F2_1 = build_T_from_D(D2, D2, W1_F2_1);
% ---------- main AC ----------
Hac_omega1_on_F1 = (abs(Om1)^2/4) * (Pa * T1_on_F1 * Pa);
Hac_omega2_on_F1 = (abs(Om2)^2/4) * (Pa * T2_on_F1 * Pa);

Hac_omega1_on_F2 = (abs(Om1)^2/4) * (Pb * T1_on_F2 * Pb);
Hac_omega2_on_F2 = (abs(Om2)^2/4) * (Pb * T2_on_F2 * Pb);

dHg_minus = Hac_omega1_on_F1 + Hac_omega2_on_F1;
dHg_plus  = Hac_omega1_on_F2 + Hac_omega2_on_F2;

% ---------- Raman ----------
H12 = (conj(Om1)*Om2/4) * exp(1i*phi) * (Pb * T12 * Pa);
Hraman = H12 + H12';

dHg_0 = Hraman;

% ---------- parts ----------
parts.AC_omega1_on_F1 = Pa * T1_on_F1 * Pa;
parts.AC_omega2_on_F1 = Pa * T2_on_F1 * Pa;
parts.AC_omega1_on_F2 = Pb * T1_on_F2 * Pb;
parts.AC_omega2_on_F2 = Pb * T2_on_F2 * Pb;

parts.Hac_omega1_on_F1 = Hac_omega1_on_F1;
parts.Hac_omega2_on_F1 = Hac_omega2_on_F1;
parts.Hac_omega1_on_F2 = Hac_omega1_on_F2;
parts.Hac_omega2_on_F2 = Hac_omega2_on_F2;
M1_raw = Pa * T1_on_F1 * Pa + Pa * T2_on_F1 * Pa;
M2_raw = Pb * T1_on_F2 * Pb + Pb * T2_on_F2 * Pb;

parts.M1 = 0.5 * (M1_raw + M1_raw');
parts.M2 = 0.5 * (M2_raw + M2_raw');
parts.M12 = Pb * T12 * Pa;

parts.H_ac_diag = dHg_minus + dHg_plus;
parts.H_raman   = Hraman;
parts.H_total   = parts.H_ac_diag + Hraman;
end

function v = getfield_or(s, name, defaultv)
    if isstruct(s) && isfield(s,name), v = s.(name); else, v = defaultv; end
end
function [W, info] = build_W_from_detune(Fp_each_excited, C, params, rule)
% Fp_each_excited : 16x1，每个激发态对应的 F'
% C              : 结构体，含 delta_e_Fp, gamma_e_Fp
% params         : 结构体，含 delta1, delta2
% rule           : 'avg' | 'delta1' | 'delta2'

    if nargin < 4 || isempty(rule)
        rule = 'avg';
    end

    dE  = C.delta_e_Fp(:).';   % 1x4
    Gam = C.gamma_e_Fp(:).';   % 1x4

    Delta1 = params.delta1 - dE;
    Delta2 = params.delta2 - dE;

    D1c = Delta1 - 1i*Gam/2;
    D2c = Delta2 - 1i*Gam/2;

    switch lower(rule)
        case 'delta1'
            wFp = 1 ./ D1c;
        case 'delta2'
            wFp = 1 ./ D2c;
        case 'avg'
            wFp = 0.5 * (1./D1c + 1./D2c);
        otherwise
            error('未知 rule: %s', rule);
    end

    W = diag(wFp(Fp_each_excited + 1));

    info = struct();
    info.Delta1 = Delta1;
    info.Delta2 = Delta2;
    info.wFp    = wFp;
end


% --------- 无权重核心：生成 DL/DR（Steck 规范，含 √2，D2 线） ---------
function [D, parts] = build_D_from_pol_auto(e_cart, quant_axis)
% 输入：
%   e_cart     : 3x1 实验室坐标偏振
%   quant_axis : 3x1 量子化轴
%
% 输出：
%   D     : 8x16，单光子耦合矩阵
%   parts : 调试信息，包括 d_q, e_q, D_parts, Fg/Fp 等

    [d_q, meta] = rb87_D2_dq_steck();

    % 复用你现有的偏振投影函数。由于这里只处理一束光，
    % 第二个输入重复传 e_cart 即可。
    pol = raman_polarizations_to_spherical(e_cart, e_cart, quant_axis);

    % 你的球基输出顺序是 [+1;0;-1]，这里重排到 [-1;0;+1]
    e_pm = pol.e_in_sph(:);
    e_q  = [e_pm(3); e_pm(2); e_pm(1)];

    [D, D_parts, dq_std] = contract_field_with_dq(d_q, e_q);

    Fg = [2*ones(5,1); 1*ones(3,1)];
    Fp = [3*ones(7,1); 2*ones(5,1); 1*ones(3,1); 0];

    parts = struct();
    parts.D       = D;
    parts.D_parts = D_parts;
    parts.dq      = dq_std;
    parts.e_q     = e_q;
    parts.pol     = pol;
    parts.Fg      = Fg;
    parts.Fp      = Fp;
    parts.meta    = meta;
end
function [D, D_parts, dq_std] = contract_field_with_dq(d_q, E_q)
% d_q 支持两种形状：
%   (3,Ng,Ne)      : 第一维是 q=-1,0,+1
%   (Ng,Ne,3)      : 第三维是 q=-1,0,+1
%
% 返回：
%   D       : Ng x Ne
%   D_parts : 3 x Ng x Ne
%   dq_std  : 统一后的 (3,Ng,Ne)

    E_q = E_q(:);
    if numel(E_q) ~= 3
        error('E_q 必须是 3x1，对应 q=-1,0,+1');
    end

    sz = size(d_q);
    if ndims(d_q) ~= 3
        error('d_q 必须是三维数组');
    end

    if sz(1) == 3
        dq_std = d_q;                  % 已经是 (3,Ng,Ne)
    elseif sz(3) == 3
        dq_std = permute(d_q, [3,1,2]);% 从 (Ng,Ne,3) -> (3,Ng,Ne)
    else
        error('d_q 的 q 维必须在第1维或第3维，且长度为3');
    end

    [~, Ng, Ne] = size(dq_std);

    D = zeros(Ng, Ne);
    D_parts = zeros(3, Ng, Ne);

    for iq = 1:3
        D_parts(iq,:,:) = conj(E_q(iq)) * squeeze(dq_std(iq,:,:));
        D = D + squeeze(D_parts(iq,:,:));
    end
end

function out = raman_polarizations_to_spherical(e_in_lab, e_ret_lab, nB)
%RAMAN_POLARIZATIONS_TO_SPHERICAL
% 将拉曼入射光和反射光的实验室偏振矢量，统一投影到
% “相对量子化轴 nB”的球基底 (sigma+, pi, sigma-)
%
% 输入:
%   e_in_lab   : 3x1 或 1x3 复偏振矢量（实验室坐标）
%   e_ret_lab  : 3x1 或 1x3 复偏振矢量（实验室坐标）
%   nB         : 3x1 或 1x3 量子化轴单位矢量（实验室坐标）
%
% 输出:
%   out 是结构体，包含
%     .e_in_sph   = [e_{+1}; e_0; e_{-1}]  入射光球分量
%     .e_ret_sph  = [e_{+1}; e_0; e_{-1}]  反射光球分量
%     .R_B        = [xB yB zB]             量子化轴局域基到实验室坐标的旋转矩阵
%     .basis      = struct('xB',xB,'yB',yB,'zB',zB)

    % 统一列向量
    e_in_lab  = e_in_lab(:);
    e_ret_lab = e_ret_lab(:);
    nB        = nB(:);

    % 归一化
    if norm(nB) < 1e-12
        error('nB 的模长太小，不能作为量子化轴。');
    end
    nB = nB / norm(nB);

    if norm(e_in_lab) < 1e-12
        error('e_in_lab 的模长太小。');
    end
    if norm(e_ret_lab) < 1e-12
        error('e_ret_lab 的模长太小。');
    end

    e_in_lab  = e_in_lab  / norm(e_in_lab);
    e_ret_lab = e_ret_lab / norm(e_ret_lab);

    %--------------------------------------------------------------
    % 构造相对量子化轴的稳定局域基 (xB, yB, zB=nB)
    % 优先用实验室 x 轴投影到垂直于 nB 的平面
    %--------------------------------------------------------------
    ex_ref = [1;0;0];
    ey_ref = [0;1;0];

    xB = ex_ref - dot(ex_ref, nB) * nB;
    if norm(xB) < 1e-10
        xB = ey_ref - dot(ey_ref, nB) * nB;
    end
    xB = xB / norm(xB);

    yB = cross(nB, xB);
    yB = yB / norm(yB);

    zB = nB;

    % 从局域基到实验室坐标的旋转矩阵
    R_B = [xB, yB, zB];

    %--------------------------------------------------------------
    % 先转到量子化轴局域笛卡尔坐标
    %--------------------------------------------------------------
    e_in_B  = R_B' * e_in_lab;
    e_ret_B = R_B' * e_ret_lab;

    %--------------------------------------------------------------
    % 局域笛卡尔 -> 球基底
    % 排序约定: [sigma+; pi; sigma-] = [q=+1; q=0; q=-1]
    %--------------------------------------------------------------
    cart_to_sph = @(v) [ ...
        -(v(1) + 1i*v(2))/sqrt(2); ...
         v(3); ...
         (v(1) - 1i*v(2))/sqrt(2) ...
    ];

    e_in_sph  = cart_to_sph(e_in_B);
    e_ret_sph = cart_to_sph(e_ret_B);

    % 输出
    out = struct();
    out.e_in_sph  = e_in_sph;
    out.e_ret_sph = e_ret_sph;
    out.R_B       = R_B;
    out.basis     = struct('xB',xB,'yB',yB,'zB',zB);
end

function T = build_T_from_D(Da, Db, W)
% Da, Db : 8x16
% W      : 16x16
% T      : 8x8
    T = Da * W * Db';
end