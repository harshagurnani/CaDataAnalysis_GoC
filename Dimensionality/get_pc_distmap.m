function [PC_Dist] = get_pc_distmap( x, ROI )

    a = ROI.Distance.pw;
    nPC = length(x.eigvec);
    PC_Dist = cell(nPC,1);
    for jj=1:nPC
        b = x.eigvec(:,jj)*x.eigvec(:,jj)'*x.eig_val(jj);
        PC_Dist{jj} = [a(:), b(:)];
    end

end