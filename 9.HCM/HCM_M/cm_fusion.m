function X = cm_fusion(LR1, LR2, HR_image, reg_param)

method  = init_option([],'method','mse_regularization');
At = reshape(LR1, [size(LR1,1)*size(LR1,2) size(LR1,3)]);
Bt = reshape(LR2, [size(LR2,1)*size(LR2,2) size(LR2,3)]);
Tt = [];
eval(['Tt = solve_T_' method '(At, Bt);']);
BB = reshape(HR_image, [size(HR_image,1)*size(HR_image,2), size(LR1,3)]);
AA=BB*Tt;
X = reshape(AA, [size(HR_image, 1), size(HR_image, 2), size(LR1, 3)]);

    function Tt = solve_T_mse(At, Bt)
        Tt = pinv(Bt' * Bt) * (Bt' * At);
    end

    function Tt = solve_T_mse_regularization(At, Bt)
        [~, s, ~] = svd(Bt, 0);
        d = diag(s);
        lamda = d(1) * reg_param;
        Tt = pinv(Bt' * Bt + lamda * eye(size(Bt, 2))) * (Bt' * At);
    end

    function val = init_option(options, name, default_val)
        if isfield(options, name), val = options.(name); else val = default_val; end
    end

end