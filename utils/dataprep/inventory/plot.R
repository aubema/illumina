catI = read.table('catI.dat')
usco = read.table('usco.dat')
uspt2 = read.table('uspt2.dat')
uswpc = read.table('uswpc.dat')

plot.dat = function(table, title) {
    plot(table[,2], table[,1], main=title, type='l')
}

pdf('cat1.pdf')
par(mfcol=c(2,2))
plot.dat(catI, 'Category I')
plot.dat(usco, 'USCO')
plot.dat(uspt2, 'USPT2')
plot.dat(uswpc, 'USWPC')
dev.off()

