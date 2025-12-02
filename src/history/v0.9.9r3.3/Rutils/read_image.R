#################################################################
#################################################################
#Snapshot (binary)
#################################################################

dir0='output' #AC2Dr output directory
lib.dir='./Rutils' #directory for lib_ac2dr.R



source(sprintf("%s/lib_ac2dr.R",lib.dir))
fl=list.files(path=dir0,pattern="_p_")

ord0=sapply(strsplit(fl, "_", fixed=TRUE), "[[", 3)
ord1=as.numeric(sapply(strsplit(ord0,".",fixed=TRUE),"[[",1))
ord2=order(ord1)

for(i in seq(from=1, by=1, to=length(ord2)))
{
    ii=ord2[i]
    fn1=sprintf('./%s/%s',dir0, fl[ii])
    aa = readbin.2dimg(fn1)
    image(aa$x, aa$y, aa$z, xlab='Angular Distance (Deg)', ylab='Elevation (m)')
    locator(1)
}

