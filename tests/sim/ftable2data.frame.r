ftable2data.frame <- function(x,...){
	y <- format(x,quote=FALSE)
	z <- data.frame(y[-1,],stringsAsFactors=FALSE)
	names(z) <- y[1,]
	z
}