function autoAdjustFigureWindowSize(figureHandle)
figureHandle.PaperPositionMode = 'Auto';
fig_pos = figureHandle.PaperPosition;
figureHandle.PaperSize = [fig_pos(3) fig_pos(4)];
end
