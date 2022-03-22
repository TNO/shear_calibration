% Assign weights to chi12 = [chi1, chi2] pairs based on weights formed in
% matrix

function w12 = assign_weight(chi12, Chi1_b, Chi2_b, W_b)

chi12_b = [Chi1_b(:), Chi2_b(:)];
w_b     = W_b(:);
nw      = size(w_b);
w12     = nan(size(chi12,1),1);

% single variale load
if all(Chi2_b == 0)
    % loop over the base (make it independent of chi12's size)
    for ii = 1:nw
        idx      = chi12(:,1) == Chi1_b(1,ii);
        
        w12(idx) = w_b(ii);
    end
% two variable loads    
else
    % loop over the base (make it independent of chi12's size)
    for ii = 1:nw
        idx1     = chi12(:,1) == chi12_b(ii,1);
        idx2     = chi12(:,2) == chi12_b(ii,2);
        idx      = idx1 & idx2;
        
        w12(idx) = w_b(ii);
    end
end

end