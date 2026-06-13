function y0_blocks = init_blocks_state(r0, v0, rho8_init)
% init_blocks_state
% 把 8x8 初始内部态密度矩阵嵌入到 40x40 block 矩阵的 block 0

Nint   = 8;
Nblock = 7;
Ntot   = Nint * Nblock;

rho40 = complex(zeros(Ntot, Ntot));

idx0 = 1:Nint;

% ---- 兼容 8x8 或 64x1 / 1x64 输入 ----
if isequal(size(rho8_init), [8,8])
    rho8 = rho8_init;
elseif numel(rho8_init) == 64
    rho8 = reshape(rho8_init(:), 8, 8);
else
    error('rho8_init must be either 8x8 or a 64-element vector.');
end

rho40(idx0, idx0) = rho8;

y0_blocks = [r0(:);
             v0(:);
             rho40(:)];
end