function [fitByX,fitByZ]=slfFindBothFits(coords,fitType)
if isempty(fitType)
    fitType='poly5';
end
Xvec=coords(1,:);Zvec=coords(3,:);
keep=~isnan(Zvec)& ~isnan(Xvec);
Zvec=Zvec(keep);Xvec=Xvec(keep);
Zvec=Zvec';Xvec=Xvec';

[fitByX,~]=fit(Xvec,Zvec,fitType);
[fitByZ,~]=fit(Zvec,Xvec,fitType);
% if statfit2.rsquare<statfit1.rsquare
%     fitFlag=1;fitA=fit1;
% else
%     fitFlag=2;fitA=fit2;
% end
end