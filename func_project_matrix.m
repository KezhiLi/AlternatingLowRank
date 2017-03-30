function [X_hat] = func_project_matrix(X, prm_struct)

if prm_struct == 1
    X_hat = func_lsqhank(X);
elseif prm_struct == 2
    X_hat = func_lsqtoep(X);
elseif prm_struct == 3
    X_hat = func_lsqcirc(X);
elseif prm_struct == 4
    X_hat = func_lsqposdef(X);
else
    disp('[LS ERROR]: unknown matrix structure')
    X_hat = X;
end

end
