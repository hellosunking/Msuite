#
# Author: Kun Sun (hellosunking@foxmail.com)
# This program is part of Msuite.
# Date: Dec 2019
#

args = commandArgs( T );
if( length(args) != 3 ) {
	print( "usage: <in.meth.log> <protocol> <out.prefix>" );
	q();
}

infile    = args[1];
outprefix = args[3];

dat = read.table( infile, head=T, comment="%" );
#chr	Reads	Total.wC	Total.wT	Total.cC	Total.cT
#chr1	75754	10064	40073	9969	40176
#chr10	44576	6061	23989	5953	23909
#chr11	45511	5882	23633	5878	23116

protocol = as.character( args[2] );
if( protocol=="TAPS" ) {
#	print( "TAPS" );
	wM = dat[,4]/(dat[,3]+dat[,4])*100;
	cM = dat[,6]/(dat[,5]+dat[,6])*100;
} else {
#	print( "BS" );
	wM = dat[,3]/(dat[,3]+dat[,4])*100;
	cM = dat[,5]/(dat[,5]+dat[,6])*100;
}
DNAm=as.matrix(rbind(wM, cM));

## PDF
out = paste( outprefix, 'pdf', sep='.' );
pdf( out, width=20, height=8 );
par( mar=c(9,5,1,0.2) );
barplot( DNAm, col=c("red", "blue"), ylim=c(0, 100), ylab="Methylation density (%)", cex.lab=1.3,
			cex.names=0.9, beside=T, xaxt="n", xlab="" );
## beside=T: group; beside=F: stack
axis(1, labels=F);
xlab.pos=(1:nrow(dat))*3;
text(xlab.pos, par("usr")[3]-2, srt=45, adj=1, labels=dat[,1], xpd=T);
legend( 'top', c("Forward chain", "Crick chain"), lty=c(1,1), col=c("red", "blue"), horiz=T, bty='n', cex=1.2 );
dev.off();

## PNG ##names.arg=dat[,1]
out = paste( outprefix, 'png', sep='.' );
png( out, width=1280, height=540 );
par( mar=c(9,5,1,0.2) );
barplot( DNAm, col=c("red", "blue"), ylim=c(0, 100), ylab="Methylation density (%)", cex.lab=1.3,
			cex.names=0.9, beside=T, xaxt="n", xlab="" );
## beside=T: group; beside=F: stack
axis(1, labels=F);
text(xlab.pos, par("usr")[3]-2, srt=45, adj=1, labels=dat[,1], xpd=T);
legend( 'top', c("Forward chain", "Crick chain"), lty=c(1,1), col=c("red", "blue"), horiz=T, bty='n', cex=1.2 );
dev.off();

