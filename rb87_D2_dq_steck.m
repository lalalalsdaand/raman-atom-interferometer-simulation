function [M_q, parts] = rb87_D2_dq_steck()
% M_q(iq, ig, ie) = <F,m | e r_q | F',m'>
% 顺序 iq: q = -1, 0, +1
% 结果按 <J=1/2||er||J'=3/2> 归一化，与 Steck 表 9–14 可直接比较

    % ground basis
    Fg = [2*ones(5,1); 1*ones(3,1)];
    mg = [ 2; 1; 0; -1; -2;  1; 0; -1];

    % excited basis
    Fp = [3*ones(7,1); 2*ones(5,1); 1*ones(3,1); 0];
    mp = [ 3; 2; 1; 0; -1; -2; -3;  2; 1; 0; -1; -2;  1; 0; -1;  0];

    Ng = numel(Fg);
    Ne = numel(Fp);

    J  = 1/2;
    Jp = 3/2;
    I  = 3/2;

    M_q = zeros(3, Ng, Ne);   % q=-1,0,+1

    for ig = 1:Ng
        F  = Fg(ig);
        m  = mg(ig);

        for ie = 1:Ne
            Fpe  = Fp(ie);
            mp_e = mp(ie);

            q = m - mp_e;   % 注意：Steck Eq.(35) 满足 m = m' + q
            if abs(q) > 1
                continue;
            end

            iq = q + 2;     % -1->1, 0->2, +1->3

            threej = Wigner3j(Fpe, 1, F, mp_e, q, -m);
            if abs(threej) < 1e-15
                continue;
            end

            red_hf = (-1)^(Fpe + J + 1 + I) ...
                     * sqrt((2*Fpe + 1)*(2*J + 1)) ...
                     * Wigner6j(J, Jp, 1, Fpe, F, I);

            M_q(iq, ig, ie) = (-1)^(Fpe - 1 + m) ...
                              * sqrt(2*F + 1) ...
                              * threej * red_hf;
        end
    end

    parts = struct('Fg',Fg,'mg',mg,'Fp',Fp,'mp',mp);
end