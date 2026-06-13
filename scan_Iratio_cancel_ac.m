function result = scan_Iratio_cancel_ac(g, C, params, idx_F1_mf0, idx_F2_mf0, ratio_list)
% 扫描 I2/I1，寻找 mF=0 两态差分 AC 光移零点
%
% idx_F1_mf0 : |F=1,mF=0> 在基底中的索引
% idx_F2_mf0 : |F=2,mF=0> 在基底中的索引

    if ~isfield(params,'I1')
        error('params.I1 必须提供，作为基准光强');
    end

    params.dipole = 3.58e-29;

    nR = numel(ratio_list);
    dAC = zeros(nR,1);
    Om1_list = zeros(nR,1);
    Om2_list = zeros(nR,1);

    c0   = 2.99792458e8;
    eps0 = 8.8541878128e-12;
    hbar = 1.054571817e-34;
    d    = params.dipole;

    for k = 1:nR
        r = ratio_list(k);

        I1 = params.Itotal/(1+r);
        I2 = params.Itotal * r/(1+r);

        params_k = params;
        params_k.I_ratio = r;

        E1 = sqrt(2*I1/(c0*eps0));
        E2 = sqrt(2*I2/(c0*eps0));

        params_k.Omega1R = d * E1 / hbar;
        params_k.Omega2R = d * E2 / hbar;

        Om1_list(k) = params_k.Omega1R;
        Om2_list(k) = params_k.Omega2R;

        [dHg_plus, dHg_minus, dHg_0,parts] = build_ac_eff_coprop(g, C, params_k);
        
        
%         preCell= init_prePair_terms(params, C);
%         parts= preCell{2,2};
       
       Hac =  (abs(Om1_list(k))^2) * parts.M1 + (abs(Om2_list(k))^2) * parts.M2 ;
%         Hac =  (abs(Om1_list(k))^2) * parts.M11 + (abs(Om2_list(k))^2) * parts.M22 ;
        % 只看 |F=2,mF=0> 与 |F=1,mF=0> 的差分 AC 光移
        dAC(k) = real(Hac(idx_F2_mf0, idx_F2_mf0) - Hac(idx_F1_mf0, idx_F1_mf0));
    end

    k0 = [];
    for k = 1:nR-1
        if dAC(k) == 0 || dAC(k)*dAC(k+1) < 0
            k0 = k;
            break;
        end
    end

    if ~isempty(k0)
        ratio_zero = interp1(dAC(k0:k0+1), ratio_list(k0:k0+1), 0, 'linear');
    else
        [~, imin] = min(abs(dAC));
        ratio_zero = ratio_list(imin);
    end

    result.ratio_list = ratio_list(:);
    result.dAC        = dAC;
    result.ratio_zero = ratio_zero;
    result.Omega1R    = Om1_list;
    result.Omega2R    = Om2_list;
    result.idx_F1_mf0 = idx_F1_mf0;
    result.idx_F2_mf0 = idx_F2_mf0;
end