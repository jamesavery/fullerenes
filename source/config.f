      module config
      integer,parameter :: Nmax=5000 !  Change NMAX if RAM is not sufficient
      integer,parameter :: msrs=56+1 !  Size of Schlegel output matrix
      integer,parameter :: NmaxL=(Nmax*(Nmax-1))/2
      integer,parameter :: Mmax=Nmax/2+2
      integer,parameter :: NSpScale=12
      integer,parameter :: NSpirals=Mmax*NSpScale
      integer,parameter :: Emax=3*Nmax/2
      integer,parameter :: maxit=2000000
      integer,parameter :: ffmaxdim=18 ! maximum number of parameters per force field
      real(8),parameter :: DPI=3.14159265358979d0 
      save
      end module config
