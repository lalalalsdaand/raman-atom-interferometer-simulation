function dydt = Raman_transition_5blocks_40(t, y, gravity, tstart, H_B, prePair, params, w_v_vec)
% Raman_transition_5blocks_40
% ------------------------------------------------------------
% 40维五块模型：
%
% block 0 : 初始动量态              |F, p>
% block X : x主通道一次跃迁后        |F, p + ħ kx>
% block Y : y主通道一次跃迁后        |F, p + ħ ky>
% block XY: 一步交叉串扰通道后       |F, p + ħ kxy>
% block C : 伴生态(先x后y/先y后x)    |F, p + ħ(kx+ky)>
%
% 每个 block 内部仍是 8 维：
%   [F=2(5态), F=1(3态)]
%
% 输入:
%   y = [r(3); v(3); rho40(:)],  rho40 为 40x40
%
% 说明:
%   这版是“五通道统一”的 40x40 截断模型，不是无限动量格点模型。
% ------------------------------------------------------------

% -------------------- 解包状态 --------------------
r = y(1:3);
v = y(4:6);

Nint   = 8;
Nblock = 7;
Ntot   = Nint * Nblock;

rho = reshape(y(7:end), Ntot, Ntot);

tp = t - tstart;

% -------------------- 常数 --------------------
hbar   = 1.054571817e-34;
mRb    = 1.44316060e-25;
lambda = 780.24e-9;
k0     = 2*pi/lambda;

% -------------------- 内部态投影 --------------------
% 基底顺序: [F=2(5), F=1(3)]
P2 = diag([1,1,1,1,1,0,0,0]);   % F=2
P1 = diag([0,0,0,0,0,1,1,1]);   % F=1

if isfield(params(1), 'omega_hf')
    omega_hf = params(1).omega_hf;   % rad/s
else
    omega_hf = 0;
end

% -------------------- 平滑开关 --------------------
tau1 = 20e-6;
A1 = 0.5*(1 + erf(tp/(sqrt(2)*tau1)));

% -------------------- block 索引 --------------------
get_idx = @(b) ((b-1)*Nint+1):(b*Nint);

idx0  = get_idx(1);   % 初始
idxX  = get_idx(7);   % x主通道
idxY  = get_idx(2);   % y主通道
idxXY = get_idx(3);   % 交叉一步通道
idxYX = get_idx(4);   % 交叉一步通道

idxX_Y  = get_idx(5);   % XY两步通道

idxY_X  = get_idx(6);   % YX两步通道
% -------------------- 缓存每个维度单光子量 --------------------
Ndim = 2;
cache = repmat(struct( ...
    'Om1',0,'Om2',0, ...
    'k1',zeros(3,1),'k2',zeros(3,1), ...
    'Delta1',0,'Delta2',0, ...
    'omega1',0,'omega2',0, ...
    'phi01',0,'phi02',0, ...
    'phi_wf_1',0,'phi_wf_2',0, ...
    'dlnG1',zeros(3,1),'dlnG2',zeros(3,1)), 1, Ndim);

c0   = 2.99792458e8;
eps0 = 8.8541878128e-12;
hbar = 1.054571817e-34;
d    =  3.58e-29;


for mu = 1:Ndim
    p = params(mu);

    % ---------- 波矢 ----------
    k1 = p.beam1.khat(:) * k0;
    k2 = p.beam2.khat(:) * (k0-1.43e2);

    % ---------- 光斑包络 ----------
    [G1, dlnG1] = gaussian_profile_xy(r, p.beam1);
    [G2, dlnG2] = gaussian_profile_xy(r, p.beam2);

    % ---------- 局域单光子 Rabi ----------
    if isfield(p,'use_I1') && p.use_I1
        I1_loc= p.Itotal /(1+p.I_ratio) ;
        I2_loc = p.Itotal * p.I_ratio /(1+p.I_ratio);
        E1 = sqrt(2*I1_loc/(c0*eps0))*G1;
        E2 = sqrt(2*I2_loc/(c0*eps0))*G2;
        Om1 =  d * E1 / hbar * A1;
        Om2 = d * E2 / hbar * A1;
    else
        Om1 = p.Omega1R * sqrt(G1) * 2*pi * A1;
        Om2 = p.Omega2R * sqrt(G2) * 2*pi * A1;
    end

    % ---------- 单光子失谐 ----------
    Delta1 = p.delta1;
    Delta2 = p.delta2;

    % ---------- 两光子扫频 ----------
    omega1 = 2*pi*w_v_vec(mu);
    omega2 = 0;

    if isfield(p, 'phi_eff')
        phi01 = p.phi_eff;
        phi02 = 0;
    else
        phi01 = 0;
        phi02 = 0;
    end

    % ---------- 波前相位 ----------
    phi_wf_1 = wavefront_phase_2d(r, tp, p.beam1, k0);
    phi_wf_2 = wavefront_phase_2d(r, tp, p.beam2, k0);

    cache(mu).Om1 = Om1;
    cache(mu).Om2 = Om2;
    cache(mu).k1  = k1;
    cache(mu).k2  = k2;
    cache(mu).Delta1 = Delta1;
    cache(mu).Delta2 = Delta2;
    cache(mu).omega1 = omega1;
    cache(mu).omega2 = omega2;
    cache(mu).phi01  = phi01;
    cache(mu).phi02  = phi02;
    cache(mu).phi_wf_1 = phi_wf_1;
    cache(mu).phi_wf_2 = phi_wf_2;
    cache(mu).dlnG1 = dlnG1;
    cache(mu).dlnG2 = dlnG2;
end

% -------------------- 有效波矢 --------------------
kx  = cache(1).k1 - cache(1).k2;   % x主通道
ky  = cache(2).k1 - cache(2).k2;   % y主通道

kab = cache(1).k1 - cache(2).k2;   % a1 - b2
kba = cache(2).k1 - cache(1).k2;   % b1 - a2

% 五块模型中 block XY 合并两条交叉路径
% 若 kab 与 kba 差异很小，这样是合理近似。
kxy = (kab );
kyx= kba;
kC = kx + ky;  % companion block 的总动量

% -------------------- 对应“格点”速度 --------------------
v0  = v;
vX  = v + hbar*kx/mRb;
vY  = v + hbar*ky/mRb;
vXY = v + hbar*kxy/mRb;
vYX = v + hbar*kyx/mRb;
vC  = v + hbar*kC/mRb;

% ============================================================
% 1) 先计算各通道（起始块出发）的一次 Raman 系数/相位/失谐
% ============================================================
[c_x0,  Phi_x0] = pair_coeff(cache(1), cache(1), v0,  omega_hf, hbar, mRb, tp );
[c_y0,  Phi_y0 ] = pair_coeff(cache(2), cache(2), v0,  omega_hf, hbar, mRb, tp );



[c_ab0, Phi_ab0] = pair_coeff(cache(1), cache(2), v0,  omega_hf, hbar, mRb, tp );
[c_ba0, Phi_ba0] = pair_coeff(cache(2), cache(1), v0,  omega_hf, hbar, mRb, tp );

% 将两条交叉路径并入一个 XY block
c_xy0    = c_ab0 ;
c_yx0    =c_ba0;
% ============================================================
% 2) companion 通道：由 X 再经 y、或由 Y 再经 x 到达
% ============================================================

[c_y_fromX, Phi_y_fromX] = pair_coeff(cache(2), cache(2), vX, omega_hf, hbar, mRb, tp );
[c_x_fromY, Phi_x_fromY] = pair_coeff(cache(1), cache(1), vY, omega_hf, hbar, mRb, tp );

% companion block 的对角失谐，取两条路径平均
[c_xy_fromY, Phi_xy_fromY] = pair_coeff(cache(1), cache(1), vXY, omega_hf, hbar, mRb, tp );

% -------------------- 光移项 --------------------
% 起始块(F=1侧)光移
a11_x = abs(cache(1).Om1)^2*exp(-1i*Phi_x0*0.00) ;
a11_y =  abs(cache(2).Om1)^2 *exp(-1i*Phi_y0*0.00);

% 一次末态(F=2侧)光移
a22_x = abs(cache(1).Om2)^2 *exp(-1i*Phi_x0*0.00) ;
a22_y = abs(cache(2).Om2)^2 *exp(-1i*Phi_y0*0.00);

%%
H22x_raw = a22_x * (P2 * prePair{1,1}.M22 * P2);
H11x_raw = a11_x * (P1 * prePair{1,1}.M11 * P1);

HX_LS=H22x_raw +H11x_raw ;
% %
% epsAC_X  = diag(HX_LS);
% deltaAC_0X = zeros(8,8);
% for ia = 1:8
%     for jb = 1:8
%         deltaAC_0X(ia,jb) = (epsAC_X(jb) - epsAC_0(ia))/hbar;
%     end
% end
%%
H22Y_raw = a22_y * (P2 * prePair{2,2}.M22 * P2);
H11Y_raw = a11_y * (P1 * prePair{2,2}.M11 * P1);
HY_LS=H22Y_raw +H11Y_raw ;

% c2x = real(trace(H22x_raw)) / 5;
% c1x = real(trace(H11x_raw)) / 3;
% 
% HY_LS =  (c2x - c1x)*(P2- P1);
%% 交叉块取两条交叉路径的平均光移

H11x_raw = a22_x * (P1 * prePair{2,1}.M11 * P1);
H22x_raw= a22_y * (P2 * prePair{1,2}.M22 * P2);

% c2x = real(trace(H22x_raw)) / 5;
% c1x = real(trace(H11x_raw)) / 3;
% 
%HXY_LS =H22x_raw +H11x_raw ;
HXY_LS=zeros(8,8);
HXY_LS(1,1)=H22x_raw(1,1);
HXY_LS(2,2)=H22x_raw(2,2)-H11x_raw(6,6);
HXY_LS(3,3)=H22x_raw(3,3)-H11x_raw(7,7);
HXY_LS(4,4)=H22x_raw(4,4)-H11x_raw(8,8);
HXY_LS(5,5)=H22x_raw(5,5);
HXY_LS(6,6)=-HXY_LS(2,2);
HXY_LS(7,7)=-HXY_LS(3,3);
HXY_LS(8,8)=-HXY_LS(8,8);
%%
H11x_raw =a22_x* (P1 * prePair{1,2}.M11 * P1) ;
H22x_raw=  a22_y* (P2 * prePair{2,1}.M22 * P2);
HYX_LS=zeros(8,8);
HYX_LS(1,1)=H22x_raw(1,1);
HYX_LS(2,2)=H22x_raw(2,2)-H11x_raw(6,6);
HYX_LS(3,3)=H22x_raw(3,3)-H11x_raw(7,7);
HYX_LS(4,4)=H22x_raw(4,4)-H11x_raw(8,8);
HYX_LS(5,5)=H22x_raw(5,5);
HYX_LS(6,6)=-HYX_LS(2,2);
HYX_LS(7,7)=-HYX_LS(3,3);
HYX_LS(8,8)=-HYX_LS(8,8);
%HYX_LS =H22x_raw +H11x_raw ;

%%
% -------------------- 构造总哈密顿量 40x40 --------------------
H = complex(zeros(Ntot, Ntot));

% ============================================================
% 对角块
% ============================================================
% block 0  : 参考块
H(idx0, idx0) = hermitize(H_B+(HX_LS+ HY_LS)/2*0);%

% block X  : 一次 x Raman 后，主要落在 F=2 流形
H(idxX, idxX) = hermitize(H_B + HX_LS*0);

% block Y  : 一次 y Raman 后，主要落在 F=2 流形
H(idxY, idxY) = hermitize(H_B + HY_LS*0);

% block XY : 一步交叉 Raman 后，主要落在 F=2 流形
H(idxXY, idxXY) = hermitize(H_B);%
% block YX : 一步交叉 Raman 后，主要落在 F=2 流形S
H(idxYX, idxYX) = hermitize(H_B);
% block C  : 二次级联后，主要回到 F=1 流形
H(idxX_Y,idxX_Y) = hermitize(H_B+ HY_LS*0);

% block C  : 二次级联后，主要回到 F=1 流形
H(idxY_X,idxY_X) = hermitize(H_B+HX_LS*0);
% ============================================================
% 非对角块 1：起始块 -> 一次块
% ============================================================
% 0 <-> X
T0X = P2 * ( c_x0 * prePair{1,1}.M12 ) * P1+P2 * ( c_y0 * prePair{2,2}.M12 ) * P1;
 H(idxX, idx0) = T0X+T0X';
H(idx0, idxX) =  (T0X+T0X')';
A=50000;
% 0 <-> Y
T0Y =P2 * ( c_ba0 * prePair{2,1}.M12 ) * P1*exp(1i*2*pi*A*t) ;
H(idxY, idxXY) = T0Y+T0Y';
H(idxXY, idxY) = (T0Y+T0Y')';
T0Y =P2 * ( c_ba0 * prePair{2,1}.M12 ) * P1*exp(-1i*2*pi*A*t) ;
H(idxY, idxYX) = T0Y+T0Y';
H(idxYX, idxY) = (T0Y+T0Y')';
% 0 <-> XY
T0XY = P2 * ( c_ab0 * prePair{1,2}.M12  ) * P1*exp(-1i*2*pi*A*t);%*exp(-1i*2*pi*50000*t)
H(idxXY, idx0) = T0XY+T0XY';
H(idx0, idxXY) =(T0XY+T0XY')';

T0XY = P2 * ( c_ab0 * prePair{1,2}.M12  ) * P1*exp(-1i*2*pi*A*t);%*exp(-1i*2*pi*50000*t)
H(idxXY, idxX) = T0XY+T0XY';
H(idxX, idxXY) = conj(T0XY+T0XY');
% 0 <-> YX
T0YX = P2 * ( c_ba0 * prePair{2,1}.M12 ) * P1*exp(1i*2*pi*A*t);%( c_ab0 * prePair{1,2}.M12  ) * P1*exp(-1i*2*pi*40000*t)+P2 * ( c_ba0 * prePair{2,1}.M12 ) * P1*exp(1i*2*pi*40000*t)
H(idxYX,idxX) = T0YX+T0YX';
H(idxX, idxYX) = (T0YX+T0YX')';

T0YX = P2 * ( c_ba0 * prePair{2,1}.M12 ) * P1*exp(1i*2*pi*A*t);%( c_ab0 * prePair{1,2}.M12  ) * P1*exp(-1i*2*pi*40000*t)+P2 * ( c_ba0 * prePair{2,1}.M12 ) * P1*exp(1i*2*pi*40000*t)
H(idxYX,idx0) = T0YX+T0YX';
H(idx0, idxYX) = (T0YX+T0YX')';


% T0YX = P2 * ( c_ba0 * prePair{2,1}.M12 ) * P1*exp(1i*2*pi*A*t)+P2 * ( c_ba0 * prePair{2,1}.M12 ) * P1*exp(-1i*2*pi*A*t);%( c_ab0 * prePair{1,2}.M12  ) * P1*exp(-1i*2*pi*40000*t)+P2 * ( c_ba0 * prePair{2,1}.M12 ) * P1*exp(1i*2*pi*40000*t)
% H(idxYX,idxXY) = T0YX+T0YX.';
% H(idxXY, idxYX) = conj(T0YX+T0YX.');
% ============================================================
% 非对角块 2：一次块 -> companion 块
% ============================================================
% X <-> C  : 再受 y 通道作用，回到 F=1
% forward主过程算符是 P2 * (...) * P1，对应 e<-g
% 此处 X->C 是反向，因此取厄米共轭%

TXC =  P2 * (c_y0  * prePair{2,2}.M12 ) * P1+P2 * (c_x0 * prePair{1,1}.M12 ) * P1;
H(idxX_Y, idxX ) = TXC+TXC';
H(idxX ,idxX_Y ) = (TXC+TXC')';


% Y <-> C  : 再受 x 通道作用，回到 F=1% +P2 * ( c_ba0 * prePair{2,1}.M12 ) * P1*exp(-1i*2*pi*0000*t)+P2 * ( c_ab0 * prePair{1,2}.M12  ) * P1*exp(1i*2*pi*0000*t)
TYC = P2 * (c_x0 * prePair{1,1}.M12 ) * P1+P2 * (c_y0  * prePair{2,2}.M12 ) * P1;
H(idxY_X, idxYX ) = TYC+TYC';
H( idxYX,idxY_X) = (TYC+TYC')';

TYC = P2 * (c_x0 * prePair{1,1}.M12 ) * P1+P2 * (c_y0  * prePair{2,2}.M12 ) * P1 ;
H(idxY_X, idxXY ) = TYC+TYC';
H( idxXY,idxY_X) =(TYC+TYC')';


% 
% TXYC =  P2 * (c_y0  * prePair{2,2}.M12 ) * P1+P2 * (c_x0 * prePair{1,1}.M12 ) * P1;
% H(idxX, idxXY ) = TXYC+TXYC.';
% H(idxXY, idxX) = conj(TXYC+TXYC.');
% % 
% TYXC = P2 * (c_y0  * prePair{2,2}.M12 ) * P1+P2 * (c_x0 * prePair{1,1}.M12 ) * P1;
% H(idxY, idxXY ) = TYXC+TYXC.';
% H(idxXY, idxY) = conj(TYXC+TYXC.');
% 
% TXYC = ( P2 * ( c_x_fromY * prePair{1,1}.M12 ) * P1 )'*0+( P2 * ( c_y_fromX  * prePair{2,2}.M12 ) * P1 )'*0;
% H(idxC, idxXY) = TXYC+TXYC';
% H(idxXY, idxC) = (TXYC+TXYC')';

% ============================================================
% 可选：若你要把 block XY 也继续耦合进 companion，需要额外指定
% 它到底是通过 x 还是 y 继续走，以及目标动量是否仍等于 kx+ky。
% 对一般几何这并不严格成立，所以这里默认不加 XY <-> C。
% ============================================================

% -------------------- 力项：先关掉，调通哈密顿量后再加 --------------------
% ------------------梯度
da11_a = a11_x * cache(1).dlnG1;
da11_b = a11_y * cache(2).dlnG1;

da22_a = a22_x * cache(1).dlnG2;
da22_b = a22_y * cache(2).dlnG2;
gradPhi_aa = zeros(3,1);
gradPhi_bb = zeros(3,1);
gradPhi_ab = zeros(3,1);

dc_aa = grad_pair_coeff(c_x0, cache(1).dlnG1, cache(1).dlnG2, gradPhi_aa);
dc_bb = grad_pair_coeff(c_y0, cache(2).dlnG1, cache(2).dlnG2, gradPhi_bb);
dc_ab = grad_pair_coeff(c_xy0, cache(1).dlnG1, cache(2).dlnG2, gradPhi_ab);

dH = complex(zeros(Ntot, Ntot, 3));

for q = 1:3
    % ---- 初态块光移梯度 ----
    dH00_LS = ...
        da11_a(q) * (P1 * prePair{1,1}.M11 * P1) + ...
        da11_b(q) * (P1 * prePair{2,2}.M11 * P1);

    % ---- 各末态块光移梯度 ----
    dHAA_LS = da22_a(q) * (P2 * prePair{1,1}.M22 * P2);
    dHBB_LS = da22_b(q) * (P2 * prePair{2,2}.M22 * P2);
    dHAB_LS = da22_b(q) * (P2 * prePair{1,2}.M22 * P2);

    % ---- Raman 非对角块梯度 ----
    dTAA = P2 * (dc_aa(q) * prePair{1,1}.M12) * P1;
    dTBB = P2 * (dc_bb(q) * prePair{2,2}.M12) * P1;
    dTAB = P2 * (dc_ab(q) * prePair{1,2}.M12) * P1;
    dTBA =P2 * (dc_ab(q) * prePair{1,2}.M12) * P1;
    % ---- 组装 dH(:,:,q) ----
    dH(idx0,  idx0,  q) = hermitize(dH00_LS);
    dH(idxX, idxX, q) = hermitize(dHAA_LS);
    dH(idxY, idxY, q) = hermitize(dHBB_LS);
    dH(idxXY, idxXY, q) = hermitize(dHAB_LS);

    dH(idxX, idx0, q) = dTAA+dTAA';
    dH(idx0,  idxX, q) =  (dTAA+dTAA')';

    dH(idxY, idx0, q) = dTBB+dTBB';
    dH(idx0,  idxY, q) = (dTBB+dTBB')';

    dH(idxXY, idx0, q) = dTAB+dTAB';
    dH(idx0,  idxXY, q) = (dTAB+dTAB')';

    dH(idxYX, idx0 , q) = dTBA+dTBA';
    dH(idx0, idxYX, q) = (dTBA+dTBA')';

    dH(idxX, idxX_Y , q) = dTAA+dTAA';
    dH(idxX_Y , idxX, q) =  (dTAA+dTAA')';

    dH(idxY, idxY_X, q) = dTBB+dTBB';
    dH(idxY_X, idxY, q) = (dTBB+dTBB')';

 
end
Fgrad = zeros(3,1);
for q = 1:3
    Fgrad(q) = -hbar * real(trace(rho * dH(:,:,q)));
end

PX  = real(trace(rho(idxX,idxX)));
PY  = real(trace(rho(idxY,idxY)));
PXY = real(trace(rho(idxXY,idxXY)));
PC  = real(trace(rho(idxY_X,idxY_X)));
P0  = real(trace(rho(idx0,idx0)));
PYX = real(trace(rho(idxYX,idxYX)));
PC2  = real(trace(rho(idxX_Y,idxX_Y)));
v_mean = P0*v0 + PX*vX + PY*vY + PXY*vXY + PC*vC+PYX*vYX+ PC2*vC;
% -------------------- 经典运动方程 --------------------
drdt = v+v_mean;
dvdt = gravity(:) + Fgrad/mRb;

% -------------------- 密度矩阵方程 --------------------
drho = -1i * (H*rho - rho*H);

% -------------------- 输出 --------------------
dydt = [drdt;
        dvdt;
        drho(:)];

end


% ============================================================
% 局部辅助函数
% ============================================================

function [cij, Phi] = pair_coeff(cache_up, cache_dn, vsite, omega_hf, hbar, mRb, t)
% pair_coeff
% ------------------------------------------------------------
% 统一计算某一对光束形成的二光子有效系数、相位、失谐
%
% cache_up : “上去”的那束光（耦合 |1> -> |e>）
% cache_dn : “下来”的那束光（耦合 |2> -> |e>）
%
% 有效拉比:
%   cij = -1/4 * Om_up^* * Om_dn * 1/2*(1/Delta_up + 1/Delta_dn) * exp(-i Phi)
%
% mode 只是为了阅读方便，不影响公式。
% ------------------------------------------------------------

kpair = cache_up.k1 - cache_dn.k2;


% 采用和你原程序一致的“theta1 - theta2”风格，
% 只是这里速度换成了 vsite
theta_up = dot(cache_up.k1, vsite) * t ...
         - cache_up.omega1* t...
         + cache_up.phi01 + cache_up.phi_wf_1;
theta_dn = dot(cache_dn.k2, vsite) * t ...
         -cache_dn.omega2 * t ...
         + cache_dn.phi02 + cache_dn.phi_wf_2;

% 可以把 tp 一并传进来。这里给出实用稳定版：
Phi = theta_up - theta_dn- hbar * dot(kpair, kpair) / (2*mRb);


% 对称分母形式
cij =  conj(cache_up.Om1) * cache_dn.Om2 ...
        * exp(-1i * Phi);

end


function Aherm = hermitize(A)
% 强制厄米化，便于数值诊断
Aherm = 0.5 * (A + A');
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
function phi_wf = wavefront_phase_2d(r, t, beam, k0)
% wavefront_phase_2d
% 返回某束光在位置 r、时刻 t 的二维波前附加相位
%
% beam 可含字段：
%   .center        [3x1]
%   .khat          [3x1]
%   .zR            瑞利长度
%   .wf_coeff_defocus
%   .wf_coeff_astig
%   .wf_coeff_xy
%   .wf_noise      标量或函数句柄

    r = r(:);
    r0 = beam.center(:);
    [ex, ey, ez] = make_beam_frame(beam.khat(:));

    d = r - r0;

    x = dot(ex, d);
    y = dot(ey, d);
    z = dot(ez, d);

    phi_curv = 0;
    if isfield(beam,'zR') && ~isempty(beam.zR) && isfinite(beam.zR)
        zR = beam.zR;
        if abs(z) < 1e-12
            Rz = inf;
        else
            Rz = z * (1 + (zR^2)/(z^2));
        end
        if isfinite(Rz)
            phi_curv = k0 * (x^2 + y^2) / (2*Rz);
        end
    end

    c_def = get_field_or_default(beam, 'wf_coeff_defocus',4e4);%4e4
    c_ast = get_field_or_default(beam, 'wf_coeff_astig',5e4);%
    c_xy  = get_field_or_default(beam, 'wf_coeff_xy',    1e3);% 

    phi_aberr = c_def*(x^2 + y^2) + c_ast*(x^2 - y^2) + c_xy*(x*y);

    phi_noise = 0;
    if isfield(beam,'wf_noise')
        if isa(beam.wf_noise,'function_handle')
            phi_noise = beam.wf_noise(t);
        else
            phi_noise = beam.wf_noise;
        end
    end

    phi_wf = phi_curv + phi_aberr + phi_noise;
end


function val = get_field_or_default(s, name, default_val)
    if isfield(s, name)
        val = s.(name);
    else
        val = default_val;
    end
end
function [ex, ey, ez] = make_beam_frame(khat)
% 为每束光建立局部正交基底 ex, ey, ez
    ez = khat(:) / norm(khat(:));

    % 选一个不平行参考轴
    if abs(dot(ez, [0;0;1])) < 0.9
        ref = [0;0;1];
    else
        ref = [1;0;0];
    end

    ex = cross(ref, ez);
    ex = ex / norm(ex);

    ey = cross(ez, ex);
    ey = ey / norm(ey);
end
function dc = grad_pair_coeff(c, dlnG_a, dlnG_b, gradPhi)
% c: 标量复耦合系数
% dlnG_a, dlnG_b, gradPhi: 3x1
%
% dc 输出 3x1，每个分量是对 x/y/z 的导数
    dc = c * (0.5*dlnG_a + 0.5*dlnG_b - 1i*gradPhi);
end