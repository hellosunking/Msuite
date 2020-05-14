#
# Author: Kun Sun (hellosunking@foxmail.com)
# This program is part of Msuite.
# Date: Dec 2019
#

args = commandArgs( T );
if( length(args) != 1 ) {
	print( "usage: <in.fqstat>" );
	q();
}

infile = args[1];
outprefix = args[1];
cols=c( "blue", "yellow", "green", "red", "black" );
#A: "azure"
#C: "Tweety bird"
#G: "Green"
#T: "Carmine"

dat = read.table( infile, head=T );
#cycle A C G T N

## normalization
for( i in 1:nrow(dat) ) {
	pass = dat[i,2] + dat[i,3] + dat[i,4] + dat[i,5];
	all  = pass + dat[i,6];

	pass = pass / 100;
	dat[i,2] = dat[i,2] / pass;
	dat[i,3] = dat[i,3] / pass;
	dat[i,4] = dat[i,4] / pass;
	dat[i,5] = dat[i,5] / pass;

	dat[i,6] = dat[i,6] / all * 100;
}
ymax=max( dat[,2],dat[,3],dat[,4],dat[,5],dat[,6] ) + 20;

out = paste( outprefix, 'pdf', sep='.' );
pdf( out, width=8, height=4 );
par( mar=c(5,5,1,1) );
plot(  dat$N ~ dat$Cycle, type='b', pch=19, col="black", ylim=c(0, ymax),
		xlab="Sequencing cycles (bp)", ylab="Frequencies (%)", cex.lab=1.5 );
lines( dat$A ~ dat$Cycle, type='b', pch=19, col=cols[1] );
lines( dat$C ~ dat$Cycle, type='b', pch=19, col=cols[2] );
lines( dat$G ~ dat$Cycle, type='b', pch=19, col=cols[3] );
lines( dat$T ~ dat$Cycle, type='b', pch=19, col=cols[4] );
legend( 'top', c('A','C','G','T','N'), col=cols, pch=rep(19, 5), horiz=T, bty='n', cex=1.5 );
dev.off();

out = paste( outprefix, 'png', sep='.' );
png( out, width=1080, height=540 );
par( mar=c(5,5,1,1) );
plot(  dat$N ~ dat$Cycle, type='b', pch=19, col="black", ylim=c(0, ymax),
		xlab="Sequencing cycles (bp)", ylab="Frequencies (%)", cex.lab=1.5 );
lines( dat$A ~ dat$Cycle, type='b', pch=19, col=cols[1] );
lines( dat$C ~ dat$Cycle, type='b', pch=19, col=cols[2] );
lines( dat$G ~ dat$Cycle, type='b', pch=19, col=cols[3] );
lines( dat$T ~ dat$Cycle, type='b', pch=19, col=cols[4] );
legend( 'top', c('A','C','G','T','N'), col=cols, pch=rep(19, 5), horiz=T, bty='n', cex=1.5 );
dev.off();
