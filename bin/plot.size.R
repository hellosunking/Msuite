#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# This program is part of Msuite.
# Date: Apr 2020
#

argv = commandArgs(T);
if( length(argv) != 1 )
{
	print( 'usage: R --slave --args <in.size> < plot.R' );
	q();
}

dat = read.table( argv[1], head=T );
#Size	Count	Proportion
#99	1741	0.194545
#100	1740	0.194433

if ( min(dat[,1]) > 100 ) {
	xmin = 100;
} else if ( min(dat[,1]) > 50 ) {
	xmin = 50;
} else {
	xmin = 0;
}

n = nrow(dat);
for ( i in 0:(n-1) ) {
	j = n - i;

	if( dat[j,3] > 0.01 ) {
		break;
	}
}

if ( dat[j,1] > 300 ) {
	xmax = dat[j,1];
} else {
	xmax = 300;
}

outfileName = paste(argv[1], "pdf", sep=".");
pdf( outfileName );
plot( dat[,3]~dat[,1], type='l', lwd=4, xlim=c(xmin, xmax),
			xlab="Fragment size (bp)", ylab="Frequency (%)", cex.lab=1.5 );
dev.off();

outfileName = paste(argv[1], "png", sep=".");
png( outfileName, width=1080, height=640 );
plot( dat[,3]~dat[,1], type='l', lwd=4, xlim=c(xmin, xmax),
			xlab="Fragment size (bp)", ylab="Frequency (%)", cex.lab=1.5 );
dev.off();


