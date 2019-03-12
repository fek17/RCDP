function figExport(w,h,name)
global printFlag
formatFig(w,h)
if printFlag == true
    print(gcf, '-dpdf', [pwd '/figures/' name '.pdf']);
end
end