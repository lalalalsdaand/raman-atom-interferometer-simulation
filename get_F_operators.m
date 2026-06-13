function F_ops = get_F_operators(g)
% 创建一个元胞数组，用三个8x8的零矩阵占位，将用于存储 Fx, Fy, Fz。
F_ops = {zeros(g.N), zeros(g.N), zeros(g.N)};


% 循环遍历8个基态中的每一个态 |r> = |F, mF>。
% 在矩阵表示中，算符作用在第 r 个基矢上，其结果体现在矩阵的第 r 列。
for r = 1:g.N
    % 提取当前基矢的量子数
    F = g.F(r); mF = g.mF(r);

    % Fz |F, mF> = mF |F, mF>  构建 Fz 矩阵
    F_ops{3}(r,r) = mF;
    % 实现F+升算符的作用
    % 升算符只能作用在 mF < F 的态上
    if mF < F
        % 计算 F+ |F, mF> 作用后的新态 |F, mF+1>
        mF_prime = mF + 1;

        % 计算系数 c = sqrt(F(F+1) - mF(mF+1))
        c = sqrt(F*(F+1) - mF*mF_prime);

        % 寻找新态 |F, mF+1> 在我们定义的8个基矢中的索引(行号)
        c_idx = find(g.F == F & g.mF == mF_prime, 1);

        if ~isempty(c_idx)
            % 根据 Fx = (F+ + F-)/2, F+ 的贡献是 c/2
            % 矩阵元 <c_idx|Fx|r> += c/2
            F_ops{1}(c_idx, r) = F_ops{1}(c_idx, r) + c/2;

            % 根据 Fy = -i/2 * (F+ - F-), F+ 的贡献是 -i*c/2
            % 矩阵元 <c_idx|Fy|r> += -i*c/2
            F_ops{2}(c_idx, r) = F_ops{2}(c_idx, r) - 1i*c/2;
        end
    end

    % 实现F-升算符的作用
    % 降算符只能作用在 mF > -F 的态上
    if mF > -F
        mF_prime = mF - 1;
        % 计算 F- |F, mF> 作用后的新态 |F, mF-1>

        c = sqrt(F*(F+1) - mF*mF_prime);
        % 计算系数 c = sqrt(F(F+1) - mF(mF-1))

        c_idx = find(g.F == F & g.mF == mF_prime, 1);
        if ~isempty(c_idx)
            % 根据 Fx = (F+ + F-)/2, F- 的贡献是 c/2
            % 矩阵元 <c_idx|Fx|r> += c/2
            F_ops{1}(c_idx, r) = F_ops{1}(c_idx, r) + c/2;
            
            F_ops{2}(c_idx, r) = F_ops{2}(c_idx, r) + 1i*c/2;
        end
    end
end
end
