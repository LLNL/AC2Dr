#!R

#R library for AC2Dr

### Read a 2d-section image
readbin.2dimg <- function (fn) {
#z is a m x n matrix, z[m][n] = z[x(1, ..., m)][y(1, ..., n)].
#x is angle in degree, y is elevation in meter.
#t is the time of image

    fptr1=file(fn, "rb")
    ftype = readBin(fptr1, integer(), n=1, size=4)
    nx1=readBin(fptr1, integer(), n=1, size=4)
    ny1=readBin(fptr1, integer(), n=1, size=4)
    dx1=readBin(fptr1, double(), n=1, size=8)
    dy1=readBin(fptr1, double(), n=1, size=8)
    bb1=matrix(readBin(fptr1, double(), n=nx1*ny1, size=8),nrow=nx1,ncol=ny1,byrow=TRUE)
    imgt = readBin(fptr1, double(), n=1, size=8)
    close(fptr1)
    xx=seq(from=0,by=dx1, length=nx1)
    yy=seq(from=0,by=dy1, length=ny1)
    aa=list(t=imgt, x=xx, y=yy, z=bb1)
    return(aa)

}


### Read a binary receiver file (waveform recorded)
readbin.rec <- function (fn) {
#p is output at time t.
#x (degree) and y (elevation) are the coordinate of receiver

    fptr=file(fn, "rb")
    ftype=readBin(fptr, integer(), n=1, size=4)
    x = readBin(fptr, double(), n=1, size=8)
    y = readBin(fptr, double(), n=1, size=8)
    dt1=readBin(fptr, double(), n=1, size=8)
    nt = readBin(fptr, integer(), n=1, size=4)
    bb=readBin(fptr, double(), n=nt, size=8)
    close(fptr)

    tt1=seq(from=0,by=dt1,length=nt)
    aa=list(x=x, y=y, t=tt1, p=bb)
    return(aa)
}


### Read a ASCII receiver file (waveform recorded)
readtxt.rec <- function (fn) {
#p is output at time t.
#x (degree) and y (elevation) are the coordinate of receiver

    aa0=scan(fn,what=0)
    x= aa0[1]
    y = aa0[2]
    dt1 = aa0[3]
    nt = aa0[4]
    bb = aa0[5:(nt+4)]
    tt1=seq(from=0,by=dt1,length=nt)
    aa=list(x=x, y=y, t=tt1, p=bb)
    return(aa)

}

### Create a 2-d material file (Binary)
makebin.2dmat <- function (x, y, z, fn, swap=0) {
#z is a m x n matrix, z[m][n] = z[x(1, ..., m)][y(1, ..., n)].
#x is angle in degree, y is elevation in meter.
   ftype=2
   dx = diff(x)[1]
   dy = diff(y)[1]
   nx = length(x)
   ny = length(y)

   fin=file(fn, "wb")
   ED=.Platform$endian
   if(swap==1) { ED="swap" }
   writeBin(as.integer(ftype), fin, size=4, endian=ED)
   writeBin(as.integer(nx), fin, size=4, endian=ED)
   writeBin(as.integer(ny), fin, size=4, endian=ED)
   writeBin(as.double(dx), fin, size=8, endian=ED)
   writeBin(as.double(dy), fin, size=8, endian=ED)

   for(i in 1:nx) {
      for(j in 1:ny) {
         writeBin(as.double(z[i,j]), fin, size=8, endian=ED)
      }
   }
   close(fin)
}

### Read a 2-d material file (Binary)
readbin.2dmat <- function (fn) {
#z is a m x n matrix, z[m][n] = z[x(1, ..., m)][y(1, ..., n)].
#x is angle in degree, y is elevation in meter.
    fptr1=file(fn, "rb")
    ftype = readBin(fptr1, integer(), n=1, size=4)
    nx1=readBin(fptr1, integer(), n=1, size=4)
    ny1=readBin(fptr1, integer(), n=1, size=4)
    dx1=readBin(fptr1, double(), n=1, size=8)
    dy1=readBin(fptr1, double(), n=1, size=8)
    bb1=matrix(readBin(fptr1, double(), n=nx1*ny1, size=8),nrow=nx1,ncol=ny1,byrow=TRUE)
    close(fptr1)
    xx=seq(from=0,by=dx1, length=nx1)
    yy=seq(from=0,by=dy1, length=ny1)
    aa=list(x=xx, y=yy, z=bb1)
    return(aa)
}

### Create a 1-d material profile (Binary)
makebin.1dprof <- function (y, z, fn, swap=0) {
#z is a vector, z[y], of elevation y (meter).
   ftype=1
   fin=file(fn, "wb")
   ED=.Platform$endian
   if(swap==1) { ED="swap" }
   writeBin(as.integer(ftype), fin, size=4, endian=ED)
   for(i in 1:length(z)) {
      writeBin(as.double(y[i]), fin, size=8, endian=ED)
      writeBin(as.double(z[i]), fin, size=8, endian=ED)
   } 
   close(fin)
}

### Read a 1-d material profile (Binary)
readbin.1dprof <- function (fn) {
#z is a vector, z[y], of elevation y (meter).
   fptr=file(fn, "rb")
   ftype = readBin(fptr, integer(), n=1, size=4)
   len=0
   while(length(readBin(fptr, double(), n=1, size=8))>0) {
      len=len+1
   }
   close(fptr)

   fptr=file(fn, "rb")
   ftype = readBin(fptr, integer(), n=1, size=4)
   aa = readBin(fptr, double(), n=len, size=8)
   close(fptr)
   ii=seq(from=1,by=2,to=length(aa))
   jj=seq(from=2,by=2,to=length(aa))
   elev=aa[ii]
   val=aa[jj]
   bb=list(y=elev, z=val)
   return(bb)
}


### Create a 1-d material profile (ASCII)
maketxt.1dprof <- function (y, z, fn) {
#z is a vector, z[y], of elevation y (meter).
   for(i in 1:length(z)) {
      if(i==1) {
         cat(sprintf("%f %f\n", y, z), file=fn, append=FALSE)
      } else {
         cat(sprintf("%f %f\n", y, z), file=fn, append=TRUE)
      }
   } 
}

### Read a 1-d material profile (ASCII)
readtxt.1dprof <- function (fn) {
#z is a vector, z[y], of elevation y (meter).
   aa=scan(fn, what=list(y=0,z=0))
   return(aa)
}




########################################################################################
### Auxiliary functions for output manipulation and display
########################################################################################
### moving average
mavg=function(x,nn) {
    y=rep(1/nn,nn)
    ny=length(y)
    nx=length(x)
    addx=floor(ny/2)
    xx=c(rep(x[1],addx),x,rep(x[nx],addx))
    res0=stats::filter(xx,y,method='convolution',sides=2)
    res=as.vector(res0)[(addx+1):(addx+nx)]
    return(res)
}

img2ras=function(A) {
    B=t(A)
    C=B[(dim(B)[1]:1),]
    return(C)
}

ras2img=function(A) {
   B=t(A)
   C=B[,dim(B)[2]:1]
   return(C)
}

#r=radius, t=theta(lat), p=phi (lon)
rtp2xyz <- function (r, t, p) {
  x = r*cos(t*pi/180)*cos(p*pi/180)
  y = r*cos(t*pi/180)*sin(p*pi/180)
  z = r*sin(t*pi/180)
  val=list(x=x, y=y, z=z)
  return(val)
}

########################################################################################
### for previous versions below v0.9.9r3
########################################################################################
### Read a 2d-section image
readbin.2dimg.v0.9.9r2 <- function (fn) {
    fptr1=file(fn, "rb")
    nx1=readBin(fptr1, integer(), n=1, size=4)
    ny1=readBin(fptr1, integer(), n=1, size=4)
    dx1=readBin(fptr1, double(), n=1, size=8)
    dy1=readBin(fptr1, double(), n=1, size=8)
    bb1=matrix(readBin(fptr1, double(), n=nx1*ny1, size=8),nrow=nx1,ncol=ny1,byrow=TRUE)
    close(fptr1)
    xx=seq(from=0,by=dx1, length=nx1)
    yy=seq(from=0,by=dy1, length=ny1)
    aa=list(x=xx, y=yy, z=bb1)
    return(aa)
}


### Read an waveform file (receiver output)
readbin.rec.v0.9.9r2 <- function (fn) {

    fptr=file(fn, "rb")
    dt1=readBin(fptr, double(), n=1, size=8)
    bb=readBin(fptr, double(), n=10^7, size=8)
    close(fptr)

    tt1=seq(from=0,by=dt1,length=length(bb))
    aa=list(t=tt1, p=bb)
    return(aa)
}

### Read a 1-d material profile (Binary)
readbin.1dprof.v0.9.9r2 <- function (fn) {
#z is a vector, z[y], of elevation y (meter).
   fptr=file(fn, "rb")
   len=0
   while(length(readBin(fptr, double(), n=1, size=8))>0) {
      len=len+1
   }
   close(fptr)

   fptr=file(fn, "rb")
   aa = readBin(fptr, double(), n=len, size=8)
   close(fptr)
   ii=seq(from=1,by=2,to=length(aa))
   jj=seq(from=2,by=2,to=length(aa))
   elev=aa[ii]
   val=aa[jj]
   bb=list(y=elev, z=val)
   return(bb)
}


