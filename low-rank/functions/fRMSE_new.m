function [RMSE, RelativeErr]=fRMSE_new(E_org, M_recover_side, Tpos)
% E_org: the real data
% M_recover_side: the predicted data
% Tpos: the position

RMSE = sqrt(sum(sum(((E_org - M_recover_side) .* Tpos).^2)) / nnz(Tpos));
RelativeErr = norm(((E_org - M_recover_side) .* Tpos), 'fro') / norm(((E_org) .* Tpos), 'fro');
end