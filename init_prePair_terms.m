function prePair = init_prePair_terms(params, C)
% init_prePair_terms
% ------------------------------------------------------------
% 生成严格二维求和所需的 prePair{mu,nu}.M11 / M22 / M12
%
% 输出:
%   prePair{mu,nu} 结构体，含
%       .M11   对应 H11^{mu,nu}
%       .M22   对应 H22^{mu,nu}
%       .M12   对应 H12^{mu,nu}
%
% 定义:
%   M11^{mu,nu} = Pa * T( beam1_mu , beam1_nu ) * Pa
%   M22^{mu,nu} = Pb * T( beam2_mu , beam2_nu ) * Pb
%   M12^{mu,nu} = Pb * T( beam1_mu , beam2_nu ) * Pa
%
% 这里 T(A,B) 表示由 A 光和 B 光构成的二阶有效算符。
%
% 说明:
%   这里调用的 dipole_T_FF_auto 已经把 detuning 放进 W 里了，
%   因此返回的是“完整二阶算符”，不是纯 CG/偏振结构。
%   如果你后面想自己再单独乘 detuning，请把 opt.use_detuning=false。
% ------------------------------------------------------------

    Ndim = numel(params);
    prePair = cell(Ndim, Ndim);

    if ~isfield(params(1), 'quant_axis')
        error('params(1).quant_axis is required.');
    end
    qa = params(1).quant_axis(:);
    qa = qa / norm(qa);

    % 8态顺序：F=2 五个态；F=1 三个态
    Fv = [2 2 2 2 2 1 1 1].';
    Pa = diag(Fv == 1);
    Pb = diag(Fv == 2);

    % 对于严格二维求和，这里统一采用 avg 规则
    % 同家族(mu,nu)会自动对应 0.5*(1/Delta_mu + 1/Delta_nu)
    % mu=nu 时自然退化成单束项
    opt = struct();
    opt.rule = 'avg';
    opt.use_detuning = true;

    for mu = 1:Ndim
        for nu = 1:Ndim

            % -------- beam1(mu) with beam1(nu): -> M11 --------
            e_mu1 = params(mu).e1_vec;
            e_nu1 = params(nu).e1_vec;
            paramsTmp = struct( ...
                'delta1', params(mu).delta1, ...
                'delta2', params(nu).delta1 );
            parts= dipole_T_FF_auto(e_mu1, e_nu1, qa, C, paramsTmp, opt);
            M11 = parts.M1;

            % -------- beam2(mu) with beam2(nu): -> M22 --------
            e_mu2 = params(mu).e2_vec;
            e_nu2 = params(nu).e2_vec;
            paramsTmp = struct( ...
                'delta1', params(mu).delta2, ...
                'delta2', params(nu).delta2 );
             parts= dipole_T_FF_auto(e_mu2, e_nu2, qa, C, paramsTmp, opt);
            M22 = parts.M2;

            % -------- beam1(mu) with beam2(nu): -> M12 --------
            paramsTmp = struct( ...
                'delta1', params(mu).delta1, ...
                'delta2', params(nu).delta2 );
            parts=dipole_T_FF_auto(e_mu1, e_nu2, qa, C, paramsTmp, opt);
            M12 =  parts.M12;

            prePair{mu,nu} = struct( ...
                'M11', M11, ...
                'M22', M22, ...
                'M12', M12 );
        end
    end
end