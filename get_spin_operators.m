function S_ops = get_spin_operators(g, I, J)

% 参考pdf 电场强度及其计算
%量子力学中的 Wigner-Eckart 定理将一个算符的矩阵元分解为两部分：
%======（1）物理部分 (约化矩阵元)：这部分与具体的磁量子数 mF 无关，
% 只和角动量的大小 F₁, F₂ 等有关。它包含了算符作用的全部“物理”信息。
%======（2）几何部分 (Wigner 3j 符号)：这部分只和角动量的“方向”
% （即磁量子数 mF）有关，它决定了跃迁的选择定则和跃迁几率的角向分布。

    % 初始化 S_ops.S, 它是一个元胞数组(cell array)，将分别存放 Sx, Sy, Sz 的 8x8 矩阵。
    S_ops.S = {zeros(g.N), zeros(g.N), zeros(g.N)};

    % 步骤一的前置准备：计算在 |J> 表象下的约化矩阵元 <J||S||J>
    % 这是一个标准结果，对于电子(S=1/2), 且基态 J=1/2, <J||S||J> = sqrt(J*(J+1)*(2*J+1))。
    reduced_me_S_on_J = sqrt(J*(J+1)*(2*J+1));


    
    % 两个 for 循环遍历 8x8 矩阵的所有元素 (row, col)。   
for row = 1:g.N            
    % row 对应末态 |F1, mF1>
    F1 = g.F(row); mF1 = g.mF(row);

    for col = 1:g.N         
        % col 对应初态 |F2, mF2>
        F2 = g.F(col); mF2 = g.mF(col);

% ============步骤一：计算 |F> 表象下的约化矩阵元 <F1||S||F2>
% 这行代码使用了角动量理论中的标准公式，通过 6j 符号进行表象变换。
% 公式: <F1||S||F2> = (-1)^(F1+J+I+1) * sqrt((2*F1+1)*(2*F2+1)) *
% {J F1 I; F2 J 1} * <J||S||J>  // { ... } 就是 Wigner 6j 符号。
        prefactor_S = (-1)^(F1 + J + I + 1) * sqrt((2*F1+1)*(2*F2+1)) * Wigner6j(J, F1, I, F2, J, 1);
        if abs(prefactor_S) < 1e-9, continue; end
        %  得到在 |F> 表象下的约化矩阵元
        reduced_me_S_on_F = prefactor_S * reduced_me_S_on_J;

% ============步骤二：应用 Wigner-Eckart 定理计算完整矩阵元
% 为了方便计算，先在球谐基矢 {S_q | q = 0, ±1} 下计算矩阵元。
% S_0 = Sz
% S_±1 分别与升降算符 S± 相关
        for q = -1:1
            % 选择定则：一个秩为 1 的张量算符 S_q 只能连接 mF 相差 q 的态。
            if mF1 ~= mF2 + q, continue; end
            
            % 计算几何因子：Wigner 3j 符号 (F1, 1, F2; -mF1, q, mF2)
            three_j = Wigner3j(F1, 1, F2, -mF1, q, mF2);
            if abs(three_j) < 1e-9, continue; end

            % Wigner-Eckart 定理的完整形式，参考电场强度等问题
            % <F1,mF1|S_q|F2,mF2> = (-1)^(F1-mF1) * 3j_symbol * <F1||S||F2>  
            me_S_q = reduced_me_S_on_F * (-1)^(F1-mF1) * three_j;

            % 最后一步：将球谐基矢下的矩阵元 (me_S_q) 转换回笛卡尔基矢 (Sx, Sy, Sz)
            if q == 0
                % S_0 就是 Sz, 直接赋值
                S_ops.S{3}(row, col) = S_ops.S{3}(row, col) + me_S_q;
            elseif q == 1
                % S_{+1} = -(Sx + iSy)/sqrt(2)
                S_ops.S{1}(row, col) = S_ops.S{1}(row, col) - me_S_q / sqrt(2);
                S_ops.S{2}(row, col) = S_ops.S{2}(row, col) + 1i * me_S_q / sqrt(2);
            elseif q == -1
                % S_{-1} = (Sx - iSy)/sqrt(2)
                S_ops.S{1}(row, col) = S_ops.S{1}(row, col) + me_S_q / sqrt(2);
                S_ops.S{2}(row, col) = S_ops.S{2}(row, col) + 1i * me_S_q / sqrt(2);
            end
        end
    end
end
end

