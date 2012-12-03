      module config
      integer,parameter :: Nmax=5000 !  Change NMAX if RAM is not sufficient
      integer,parameter :: msrs=56+1 !  Size of Schlegel output matrix
      integer,parameter :: NmaxL=(Nmax*(Nmax-1))/2
      integer,parameter :: Mmax=Nmax/2+2
      integer,parameter :: NSpScale=12
      integer,parameter :: NSpirals=Mmax*NSpScale
      integer,parameter :: Emax=3*Nmax/2
      integer,parameter :: maxit=2000000
      integer,parameter :: maxRS=20000000
      integer,parameter :: LimitAll=122
      integer,parameter :: LimitIPR=170
      integer,parameter :: ffmaxdim=19 ! maximum number of parameters per force field
      integer,parameter :: intmax8=2147483647 ! maximum integer before overflow
      real(8),parameter :: dpi=3.14159265358979d0 
      real(8),parameter :: deg2rad=dpi/1.8d2
      real(8),parameter :: rad2deg=1.8d2/dpi
      real(8),parameter :: au2wavenumbers=2.19474625d5
      real(8),parameter :: au2eV=2.72113957d1
      real(8),parameter :: dynpercm2auperaa=2.29371049d-6
      real(8),parameter :: au2angstroem=.529177249d0
      save
      end module config
