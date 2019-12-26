function fitLine=slfCreateFitline(fitA,fitFlag,coords,hemi)
p1=fitA.p1;p2=fitA.p2;p3=fitA.p3;p4=fitA.p4;p5=fitA.p5;p6=fitA.p6;
if fitFlag==1
    fitLine=linspace(min(coords(1,:)),max(coords(1,:)),500);
    if strcmp(hemi,'left')
        fitLine=flip(fitLine,2);
    end
elseif fitFlag==2
    fitLine=linspace(min(coords(3,:)),max(coords(3,:)),500);
    fitLine=flip(fitLine,2);
end

fitLine(2,:)= p1*fitLine.^5 + p2*fitLine.^4 + p3*fitLine.^3 + p4*fitLine.^2 + p5*fitLine + p6;
if fitFlag==2
    fitLine=flip(fitLine);
end

end