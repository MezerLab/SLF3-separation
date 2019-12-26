function dice = slfGetDiceStreamline(fg,slf3IndicesParam)
% slfGetDiceStreamline reutrns the Dice coefficient calculated based on
% streamline-wise agreement between SLF3 as identified by TractSeg (fg(3))
% and SLF3 as identified by some local parameter (e.g., FA or T1)

origClass=[ones(length(fg(1).fibers),1);2*ones(length(fg(2).fibers),1);3*ones(length(fg(3).fibers),1)];
    
slf3IndicesTractseg = find(origClass==3);
overlap = length(find(ismember(slf3IndicesParam,slf3IndicesTractseg)));
paramFiberNum = length(slf3IndicesParam);
tractsegFiberNum = length(slf3IndicesTractseg);
dice = 2*overlap/(tractsegFiberNum + paramFiberNum);