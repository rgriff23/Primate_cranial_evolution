plotWireframe = function (A, links, params=NULL) 
{
  open3d()
  if (!is.null(params)) {par3d(params)}
  rgl.bg(sphere=TRUE, color=c("black"), lit=FALSE, back="fill")
  A3d = matrix(A, dim(A)[[1]], dim(A)[[2]], dimnames=dimnames(A)[1:2])
  plot3d(A3d, type = "s", col = "white", xlab = "", ylab = "", zlab = "", size = 1, aspect = FALSE, box=FALSE, axes=FALSE)
  points3d(A, color = "white", size = 4)
  for (i in 1:nrow(links)) {
    segments3d(rbind(A[links[i,1],,], A[links[i,2],,]), lwd = 2, col="green")
  }
}

