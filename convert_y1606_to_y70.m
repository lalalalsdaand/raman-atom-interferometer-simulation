function y_small = convert_y1606_to_y70(y_big)
% convert_y1606_to_y70
% ------------------------------------------------------------
% 输入:
%   y_big : Nt x 1606
%           第1~3列  : r
%           第4~6列  : v
%           第7~1606 : rho40(:)，其中 rho40 是 40x40
%
% 输出:
%   y_small : Nt x 70
%             第1~3列 : r
%             第4~6列 : v
%             第7~70  : rho8(:)，其中 rho8 是 8x8 约化密度矩阵
% ------------------------------------------------------------

    Nt = size(y_big,1);

    Nint   = 8;
    Nblock = 7;
    Ntot   = Nint * Nblock;   % 40

    if 6 + Ntot*Ntot ~= (Nint*Nblock)^2+6  
        error('Dimension mismatch: expected 6 + 40*40 = 1606.');
    end

    % 必须用复数数组
    y_small = complex(zeros(Nt, 70));

    % 位置和速度直接拷贝
    y_small(:,1:6) = y_big(:,1:6);

    % block 索引
    idx0  = 1:8;
    idxAA = 9:16;
    idxBB = 17:24;
    idxAB = 25:32;
    idxBA = 33:40;
    idxBC = 41:48;
    idxCB = 49:56;
    for it = 1:Nt
        % 按列向量口径还原 rho40
        rho40 = reshape(y_big(it,7:end).', Ntot, Ntot);

        % 对动量块做偏迹，得到 8x8 内部态约化密度矩阵
        rho8 = rho40(idx0,  idx0 ) ...
             + rho40(idxAA, idxAA) ...
             + rho40(idxBB, idxBB) ...
             + rho40(idxAB, idxAB) ...
             + rho40(idxBA, idxBA)...
             + rho40(idxBC, idxBC)...
             + rho40(idxCB , idxCB );

        % 按列展开成 64 个元素
        y_small(it,7:end) = rho8(:).';
    end
end