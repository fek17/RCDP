function figExport(w,h,name)
global t
formatFig(w,h)
if t.export == true
    print(gcf, '-dpdf', [pwd '/figures/' name '.pdf']);
end
end