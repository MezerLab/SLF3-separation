cmap = [141,211,199; 
253,192,134;%255,255,179;
190,186,218;
251,128,114;
128,177,211;
253,180,98;
179,222,105;
144 160 80]./255;

% 1 = FA
% 2 = T1
% 3 = MTV
% 4 = MD
% 5 = T2* (T2 is cmap(5,:)*0.7)
% 8 = T2w/T1w

%% TractSeg colors for SLF 1,2,3, based on Catani's colormap
cmapTs = [138,198,183;36,32,192;143,61,181]/255;
cmapTs(2,:) = cmapTs(2,:)*1.2;
cmapTs(3,:) = cmapTs(3,:)*1.2;

round(cmapTs*255)
