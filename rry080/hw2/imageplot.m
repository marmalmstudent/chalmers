function imageplot(input)

colormap(gray)
imagesc(sqrt(sqrt(abs(input))))
brighten(.2)

