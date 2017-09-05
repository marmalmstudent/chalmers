function imageplot(input)

colormap(gray)
imagesc(sqrt(sqrt(abs(input))))
% brighten(.2)
axis([0, size(input, 2)-1, 0, size(input, 1)-1])