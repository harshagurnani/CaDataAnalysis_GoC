function slope = get_v_slope( cc, ROI, vec )

    nROI = size(cc, 1);
    dist_v = [];
    cc_v = [];
    for roi1 = 1:nROI
        for roi2=roi1+1:nROI
            rvec = [(ROI.Coordinates.ROIs_X(roi1,1)-ROI.Coordinates.ROIs_X(roi2,1)),  ...
                    (ROI.Coordinates.ROIs_Y(roi1,1)-ROI.Coordinates.ROIs_Y(roi2,1))];
            dist_v = [dist_v; abs( dot(rvec, vec/norm(vec)) )];
            cc_v = [cc_v; cc(roi1, roi2)];
        end
    end
    scatter(dist_v, cc_v);hold on;
    slope = dist_v\ cc_v;

end