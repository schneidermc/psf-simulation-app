
function [stack_cent, m_focus] = centering(stack)

    stack_norm = stack ./ sum(sum(stack,1),2); 
    [Nx, Ny, ~] = size(stack);
    [~, maxidx] = max(stack_norm(:));
    [maxrow, maxcol, m_focus] = ind2sub(size(stack), maxidx);
 
    %centering the PSF within the frame
    stack_cent = circshift(stack, round(([Nx Ny]+1)/2)-[maxrow, maxcol]); 