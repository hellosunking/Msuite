#
# Author: Kun Sun (hellosunking@foxmail.com)
# This program is part of Msuite.
# Date: Dec 2019
#

args = commandArgs( T );
if( length(args) != 2 ) {
	print( "usage: <in.mbias> <mode=BS|TAPS>" );
	q();
}

mbias = read.table( args[1], head=T, comment="%" );

## calculate DNAm
if( args[2] == "BS" ) {
	wm = mbias$wC/(mbias$wC+mbias$wT) * 100;
	cm = mbias$cC/(mbias$cC+mbias$cT) * 100;
} else {
	wm = mbias$wT/(mbias$wC+mbias$wT) * 100;
	cm = mbias$cT/(mbias$cC+mbias$cT) * 100;
}

out = paste( args[1], 'pdf', sep='.' );
pdf( out, width=8, height=4 );
par( mar=c(5,5,1,1) );
plot(  wm ~ mbias$Cycle, type='o', pch=1, col="red", ylim=c(0, 100),
		xlab="Sequencing cycles (bp)", ylab="Methylation level (%)", cex.lab=1.5 );
lines( cm ~ mbias$Cycle, type='o', pch=8, col="blue" );
legend( 'top', c('Watson', 'Crick'), col=c("red", "blue"), pch=c(1, 8), horiz=T, bty='n', cex=1.5 );
dev.off();

out = paste( args[1], 'png', sep='.' );
png( out, width=1080, height=540 );
par( mar=c(5,5,1,1) );
plot(  wm ~ mbias$Cycle, type='o', pch=1, col="red", ylim=c(0, 100),
		xlab="Sequencing cycles (bp)", ylab="Methylation level (%)", cex.lab=1.5 );
lines( cm ~ mbias$Cycle, type='o', pch=8, col="blue" );
legend( 'top', c('Watson', 'Crick'), col=c("red", "blue"), pch=c(1, 8), horiz=T, bty='n', cex=1.5 );
dev.off();

