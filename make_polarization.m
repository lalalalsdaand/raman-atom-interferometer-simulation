function [e_hat, jones, basis, k_hat] = make_polarization(k_init, rot_angles_deg, orientation_deg, ellipticity_deg)
% MAKE_POLARIZATION - 根据初始波矢和多轴旋转 + 偏振参数生成偏振矢量
%
% 输入:
%   k_init           : 初始波矢 (3×1)，通常取 [1;0;0] 表示沿 +x
%   rot_angles_deg   : [rx, ry, rz] 三个旋转角度 (度)，依次绕 x,y,z 轴旋转
%   orientation_deg  : 偏振椭圆长轴方向角 ψ (度)
%   ellipticity_deg  : 偏振椭圆角 χ (度)，0=线偏振，±45=圆偏振
%
% 输出:
%   e_hat : 3×1 复数偏振单位矢量
%   jones : 2×1 复数 Jones 向量
%   basis : 3×2 横向实数单位正交基
%   k_hat : 3×1 传播方向

    % -------- 初始波矢 --------
    k_hat = k_init(:)/norm(k_init);

    % -------- 旋转矩阵 (欧拉角，Z-Y-X 顺序或你需要的顺序) --------
    rx = deg2rad(rot_angles_deg(1));
    ry = deg2rad(rot_angles_deg(2));
    rz = deg2rad(rot_angles_deg(3));

    Rx = [1 0 0;
          0 cos(rx) -sin(rx);
          0 sin(rx) cos(rx)];
    Ry = [cos(ry) 0 sin(ry);
          0 1 0;
          -sin(ry) 0 cos(ry)];
    Rz = [cos(rz) -sin(rz) 0;
          sin(rz)  cos(rz) 0;
          0        0       1];

    % 总旋转：Rz * Ry * Rx (可换顺序，看实验几何定义)
    R = Rz * Ry * Rx;

    % 旋转后的波矢
    k_hat = R * k_hat;
    k_hat = k_hat / norm(k_hat);

    % -------- 构造横向基 --------
    ref = [0;0;1];
    if abs(dot(ref,k_hat)) > 0.98
        ref = [0;1;0];
    end
    u1 = cross(k_hat, ref); u1 = u1 / norm(u1);
    u2 = cross(k_hat, u1);  u2 = u2 / norm(u2);
    basis = [u1 u2];

    % -------- Jones 矢量 --------
    psi = deg2rad(orientation_deg);
    chi = deg2rad(ellipticity_deg);
    j1 = cos(psi)*cos(chi) - 1i*sin(psi)*sin(chi);
    j2 = sin(psi)*cos(chi) + 1i*cos(psi)*sin(chi);
    jones = [j1; j2];
    jones = jones / norm(jones);

    % -------- 偏振向量 --------
    e_hat = basis * jones;
    e_hat = e_hat / norm(e_hat);
end
