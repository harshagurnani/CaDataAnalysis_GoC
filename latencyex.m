for jj=1:112
[pks] = find( zscored_puff_dff(jj,14:40) >= 1.96,1);
if ~isempty(pks)
[pks2] = find( zscored_puff_dff(jj,13+pks:70) < 1.5);
if ~isempty(pks2)
pkid6(jj,1:3)=[pks+13, pks2(1)+pks+13, pks2(end)+pks+13];
else
pkid6(jj,1:3)=[pks+13 70 70];
end
else
pkid6(jj,1:3)=[NaN NaN NaN];
end
end

early_start = pkid6(:,1);
last_respindx = pkid6(:,2);
for jj=1:112
if ~isnan(early_start(jj))
early_start(jj) = pufftimes( early_start(jj));
end
if ~isnan(last_respindx(jj))
last_respindx(jj) = pufftimes(last_respindx(jj));
end
end
resp_range = (last_respindx - early_start)
rise_times=pk_indx-early_start