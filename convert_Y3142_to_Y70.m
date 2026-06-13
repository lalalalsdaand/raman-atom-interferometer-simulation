function Y70 = convert_Y3142_to_Y70(Y56, mode, block_id, do_normalize)
% convert_Y3142_to_Y70
% ------------------------------------------------------------
% 将 Y56: [Ns, 3142, Nt]
% 转为 Y70: [Ns, 70, Nt]
%
% 其中：
%   Y56(:,1:3,:)   = 速度
%   Y56(:,4:6,:)   = 位置
%   Y56(:,7:end,:) = rho56(:), 56*56=3136
%
% 输出：
%   Y70(:,1:3,:)   = 速度
%   Y70(:,4:6,:)   = 位置
%   Y70(:,7:end,:) = rho8(:), 8*8=64
%
% mode:
%   'sum'   : 对 7 个 8x8 block 求和（推荐，物理上是对动量偏迹）
%   'block' : 取指定 block_id 的 8x8 block
%   'first' : 取第一个 8x8 block
%
% block_id:
%   当 mode='block' 时需要，范围 1~7
%
% do_normalize:
%   是否对 rho8 归一化，默认 false
% ------------------------------------------------------------

    if nargin < 2 || isempty(mode)
        mode = 'sum';
    end

    if nargin < 3 || isempty(block_id)
        block_id = [];
    end

    if nargin < 4 || isempty(do_normalize)
        do_normalize = false;
    end

    s = size(Y56);
    if numel(s) ~= 3
        error('输入 Y56 必须是 [Ns, 3142, Nt] 三维数组');
    end

    Ns = s(1);
    Ny = s(2);
    Nt = s(3);

    if Ny ~= 3142
        error('第二维必须为 3142 (= 6 + 56*56)');
    end

    Y70 = zeros(Ns, 70, Nt, 'like', Y56);

    % 保留速度和位置
    Y70(:,1:6,:) = Y56(:,1:6,:);

    for is = 1:Ns
        for it = 1:Nt

            % 取出 rho56 向量并还原成 56x56
            rho56_vec = reshape(Y56(is, 7:end, it), [], 1);
            rho56 = reshape(rho56_vec, 56, 56);

            % 压缩成 rho8
            rho8 = rho56_to_rho8_local(rho56, mode, block_id, do_normalize);

            % 写回
            Y70(is, 7:end, it) = reshape(rho8, 1, 64);
        end
    end
end


function rho8 = rho56_to_rho8_local(rho56, mode, block_id, do_normalize)
% 局部子函数：56x56 -> 8x8

    if nargin < 2 || isempty(mode)
        mode = 'sum';
    end

    if nargin < 3 || isempty(block_id)
        block_id = [];
    end

    if nargin < 4 || isempty(do_normalize)
        do_normalize = false;
    end

    rho8 = zeros(8,8, 'like', rho56);

    switch lower(mode)

        case 'first'
            rho8 = rho56(1:8,1:8);

        case 'block'
            if isempty(block_id)
                error('mode=block 时必须提供 block_id');
            end
            if block_id < 1 || block_id > 7
                error('block_id 必须在 1~7');
            end
            idx = (block_id-1)*8 + (1:8);
            rho8 = rho56(idx,idx);

        case 'sum'
            for b = 1:7
                idx = (b-1)*8 + (1:8);
                rho8 = rho8 + rho56(idx,idx);
            end

        otherwise
            error('未知 mode');
    end

    % 数值上强制厄米
    rho8 = 0.5 * (rho8 + rho8');

    % 可选归一化
    if do_normalize
        tr = trace(rho8);
        if abs(tr) > 0
            rho8 = rho8 / tr;
        end
    end
end