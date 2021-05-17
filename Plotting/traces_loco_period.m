figure;
pp2=get_loco_period(temp,true,305,[0.5,0.05],400);
pp=pp2;
min_pp=min(pp(:,3));
max_pp=max(pp(:,3));
r_pp=max_pp-min_pp
figure;
for jj=1:size(pp,1)
hold on
a=area(pp(jj,1:2),repmat(10,[1,2]));
a.FaceAlpha=0.3;
a.EdgeAlpha=0;
a.FaceColor= [ (pp(jj,3)-min_pp)/r_pp, 0 , 1-(pp(jj,3)-min_pp)/r_pp];
end
hold on
for jj=1:12
plot(spd{jj}(:,1),-spd{jj}(:,2)/5+7, 'Color', cmap(mod(jj,7)+1,:))
end
ylim=[0 10]
plot_traces(-n*10,t,[1 3 2])