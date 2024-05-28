#################################################################
#################################################################
#Waveform, pressure (binary)
#################################################################

fn='./output/STA05.bin' #AC2Dr receiver waveform file
lib.dir='./Rutils' #directory for lib_ac2dr.R


source(sprintf("%s/lib_ac2dr.R",lib.dir))
aa=readbin.rec(fn)
plot(aa$t, aa$p, type='l', xlab='Time (s)', ylab='Pressure (Pa)')


