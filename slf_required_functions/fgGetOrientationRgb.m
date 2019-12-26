function rgb = fgGetOrientationRgb(fg)
% fgGetOrientationRgb returns a cell array with the RGB colors of each node
% of each streamline.
%
% E.g.: AFQ_RenderFibers(fgTmp,'color',fgGetOrientationRgb(fg),'tube',1,'radius',[0.2 0.3])

rgb = cellfun(@(x) bsxfun(@rdivide, abs([gradient(x(1,:)'), gradient(x(2,:)'), gradient(x(3,:)')]'), vecnorm([gradient(x(1,:)'), gradient(x(2,:)'), gradient(x(3,:)')]'))', fg.fibers, 'UniformOutput', false);
