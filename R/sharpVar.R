#' Calculates variance according to Arronow...
#' @return The function returns XXX
#' @export 
#' 
#' 
#' 

#### A sharp variance estimate 
# Inputs:
#	yt: observed outcome under treatment
#	yc: observed outcome under control
# Outputs:
#	variance estimate

sharpVar<-function(yt,yc,N=length(c(yt,yc)),upper=TRUE){
	m<-length(yt)
	n<-m+length(yc)
	FPvar<-function(x,N) (N-1)/(N*(length(x)-1)) * sum((x-mean(x))^2)
	yt<-sort(yt)
	if(upper==TRUE) yc<-sort(yc) else 
		yc<-sort(yc,decreasing=TRUE)
	p_i<-unique(sort(c(seq(0,n-m,1)/(n-m),seq(0,m,1)/m))) - .Machine$double.eps^.5
	p_i[1]<-.Machine$double.eps^.5
	yti<-yt[ceiling(p_i*m)]
	yci<-yc[ceiling(p_i*(n-m))]
	p_i_minus<-c(NA,p_i[1:(length(p_i)-1)])
	return(((N-m)/m * FPvar(yt,N) + (N-(n-m))/(n-m) * FPvar(yc,N) 
		+ 2*sum(((p_i-p_i_minus)*yti*yci)[2:length(p_i)])
		- 2*mean(yt)*mean(yc))/(N-1))	
}









