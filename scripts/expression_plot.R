args <- commandArgs(TRUE)

expression_data <- read.table(args[1], header=T);

position = c(head(expression_data$position,n=1),expression_data$position,tail(expression_data$position,n=1))
expression = c(0,expression_data$expression,0);

pdf(args[2]);

limx <- c(min(position),max(position));
limy <- c(min(expression),max(expression));

plot(position,expression,type="n",xlim=limx,ylim=limy,xlab="Position",ylab="Expression");
lines(position,expression);

if (length(args) > 3)
{
	breakpos <- as.integer(args[3]);
	strand <- as.integer(args[4]);
	
	segments(breakpos, limy[1], breakpos, limy[2], col="red", lwd=3);
	
	mid <- (limy[2] + limy[1]) / 2;
	arrowlen = (limx[2] - limx[1]) * 0.1;

	arrowstart <- breakpos - arrowlen;
	if (strand < 0)
	{
		arrowstart <- breakpos + arrowlen;
	}
	
	arrows(arrowstart, mid, breakpos, mid, col="red", lwd=3);
}

dev.off();
