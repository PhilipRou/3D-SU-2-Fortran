program su2_3D

!!!                                                                                   {this option below is mine}
!!! gfortran -cpp -std=f2018 -fimplicit-none -ffree-line-length-none -fall-intrinsics -fmax-stack-var-size=65535 -fmax-errors=1 -Wall -Wextra -Warray-temporaries -fopenmp  -Ofast -o su2_3D su2_3D.f90 # -fprotect-parens -static-libgfortran -fcheck=all
!!! gfortran -cpp -std=f2018 -fimplicit-none -ffree-line-length-none -fall-intrinsics -fmax-stack-var-size=65535 -fmax-errors=1 -Wall -Wextra -Warray-temporaries -fopenacc -Ofast -o su2_3D su2_3D.f90 # -fprotect-parens -static-libgfortran -fcheck=all

!!! module load intel
!!! ifort -fpp -stand f08 -implicitnone -diag-disable 5268 -warn all -nogen-interfaces -O2 -xcore-avx2   -qopenmp -o su2_3D su2_3D.f90 # -static -static-intel
!!! ifort -fpp -stand f08 -implicitnone -diag-disable 5268 -warn all -nogen-interfaces -O2 -xcore-avx512 -qopenmp -o su2_3D su2_3D.f90

!!! module load NAGWare_f95
!!! nagfor -fpp -D__NAG__ -f2018 -free -maxcontin=1024 -u -thread_safe -kind=byte -O3 -target=core2 -openmp -o su2_3D su2_3D.f90 # -unsharedrts

!!! module load pgi
!!! module load NVHPC
!!! pgfortran -cpp -D __PGI__ -mp                   -Minfo=par   -O3 -o su2_3D su2_3D.f90 # -Bstatic -tp=skylake -Mcache_align -Mpre -Mdse -Mfma -Minfo=all
!!! pgfortran -cpp -D __PGI__ -acc -ta=nvidia,ccall -Minfo=accel -O3 -o su2_3D su2_3D.f90 # -Bstatic -tp=skylake -ta=nvidia,ccall,fastmath,fma,deepcopy

!!! export OMP_PLACES=threads OMP_PROC_BIND=true OMP_DYNAMIC=false; unset OMP_PLACES OMP_PROC_BIND OMP_DYNAMIC
!!! export OMP_PLACES=cores   OMP_PROC_BIND=true OMP_DYNAMIC=false; unset OMP_PLACES OMP_PROC_BIND OMP_DYNAMIC
!!! export OMP_PLACES=cpus    OMP_PROC_BIND=true OMP_DYNAMIC=false; unset OMP_PLACES OMP_PROC_BIND OMP_DYNAMIC

#ifdef _OPENMP
   use omp_lib, only: omp_get_num_procs,omp_get_max_threads,omp_set_num_threads
#endif

   implicit none
   integer,parameter :: sp=kind(0.0E0)
   integer,parameter :: dp=kind(0.0D0)
   integer,parameter :: Nx=32 ,Ny=Nx, Nz=Ny
   real(kind=sp)     :: beta_vals(7)
   integer(kind=sp)  :: beta_ind
   !real(kind=sp),parameter :: twopi_sp=8.0_sp*atan(1.0_sp)
   !real(kind=dp),parameter :: twopi_dp=8.0_dp*atan(1.0_dp)

   beta_vals=[6.00,7.50,9.00,10.50,12.00,13.50,15.00]
   do beta_ind=1,size(beta_vals)
      call su2_3D_main(beta_vals(beta_ind))
   end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine su2_3D_main(beta_3D)

      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz) :: U,V
      real(kind=sp), intent(in) :: beta_3D
      real(kind=sp) :: hitsize, tmp_metro, tmp_ovlax, accrate_metro, accrate_ovlax
      real(kind=dp) :: swil_thin, swil_smth, sopt_thin, sopt_smth
      real(kind=sp) :: t_sec
      integer(kind=sp) :: n_ther, n_meas, n_sepa, n_over, n_hit, n_stout, loop_ther, loop_meas, loop_sepa, loop_over, n_procs, n_threads
      integer(kind=dp) :: t_ini, t_fin, t_old, t_new, t_min, clock_rate
      character(len=80) :: hostname, hostdate, hosttime, hostzone, str_places, str_procbind, str_dynamic, str_kmp
      integer(kind=sp), allocatable :: smearlist(:)
      integer(kind=sp) :: z_loop, write_loop
      real(kind=sp) :: rho_2D , rho_3D !, rho_4D !!! rho_2D<0.24, rho_3D<0.16, rho_4D<0.12
      real(kind=dp) :: swil_corr_thin(Nz), swil_corr_smth(Nz)
      real(kind=dp), allocatable :: swil_smeared(:), swil_cross_corr(:,:,:)
      character(len=80) :: beta_str, Nz_str, Nx_str, n_stout_str, rho_str !, sim_count
      character(len=512) :: base_name, log_file_name, timeseries_file_name, corr_file_name, swil_smeared_file_name, crosscorr_file_name
#ifdef _OPENMP
      logical,parameter :: do_openmp=.true.
#else
      logical,parameter :: do_openmp=.false.
#endif

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! call su2_3D_init(U)
      ! do loop_ther = 1,200
      !    call su2_3D_metropolis  (U,6.0,0.1,4,tmp_metro)
      ! end do
      
      ! call su2_2D_stout(U,V,0.24)
      ! print *, "U(:,3,1,1,1) = ", U(:,3,1,1,1)
      ! print *, "V(:,3,1,1,1) = ", V(:,3,1,1,1)
      ! stop
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (mod(Nx,2)+mod(Ny,2)+mod(Nz,2)>0) then
         write(*,*) 'su2_3D_main: each one of Nx,Ny,Nz must be even'; stop
      end if

      ! beta_3D=8.0
      hitsize=0.1
      rho_3D =0.16 
      rho_2D =0.24 

      n_ther=100
      n_meas=5000
      n_sepa= 10
      n_over=  4
      n_hit =  8
      n_stout= 7
      allocate(smearlist(1:5))
      smearlist = [0,1,2,5,7]

      write(beta_str, "(F5.2)") beta_3D; beta_str = trim(adjustl(beta_str))
      write(Nz_str, "(I3)") Nz; Nz_str = trim(adjustl(Nz_str))
      write(Nx_str, "(I3)") Nx; Nx_str = trim(adjustl(Nx_str))
      write(n_stout_str, "(I2)") n_stout; n_stout_str = trim(adjustl(n_stout_str))
      write(rho_str, "(F5.2)") rho_2D; rho_str = trim(adjustl(rho_str))
      base_name = adjustl("beta_"//trim(beta_str)//"_Nz_"//trim(Nz_str)//"_Nx_"//trim(Nx_str)//"_n_stout_"//trim(n_stout_str)//"_rho_2D_"//trim(rho_str))
      if (allocated(smearlist)) then
         base_name = adjustl("beta_"//trim(beta_str)//"_Nz_"//trim(Nz_str)//"_Nx_"//trim(Nx_str)//"_smearlist_rho_2D_"//trim(rho_str))
      endif
      timeseries_file_name    = trim(base_name)//".txt"
      log_file_name           = trim(base_name)//"_log.txt"
      corr_file_name          = trim(base_name)//"_corr.txt"
      crosscorr_file_name     = trim(base_name)//"_crosscorr.txt"
      swil_smeared_file_name  = trim(base_name)//"_swil_smeared.txt"
      
      open(10, file = log_file_name, status = "new")
      write(10, "(a)") "Simulation parameters:"
      write(10, "(a)") " "
      write(10, "(a)") "Comment: beta scan to find parameters with good coupling to glueball correlator"
      write(10, "(a, F7.4)") "beta_3D   = ", beta_3D
      write(10, "(a, F7.4)") "rho_2D    = ", rho_2D
      write(10, "(a, F7.4)") "rho_3D    = ", rho_3D
      write(10, "(a, I7)")   "n_ther    = ", n_ther
      write(10, "(a, I7)")   "n_meas    = ", n_meas
      write(10, "(a, I7)")   "n_sepa    = ", n_sepa
      write(10, "(a, I7)")   "n_over    = ", n_over
      write(10, "(a, I7)")   "n_hit     = ", n_hit
      if (allocated(smearlist)) then
         write(10, "(a)", advance="no") "smearlist = ["
            do write_loop = 1,size(smearlist)
               write(10, "(I2,a)", advance="no") smearlist(write_loop), ", "
            end do
         write(10, "(a)") "]"
      else 
         write(10, "(a, I7)")   "n_stout   = ", n_stout
      endif
      write(10, "(a, F7.4)") "hitsize prior therm. = ", hitsize
      close(10)
      
      open(11, file = timeseries_file_name, status = "new")
      write(11, "(a)") "s_wil    s_opt    smeared_s_wil  smeared_clover  accrate_metro  accrate_ovlax  time"
      close(11)

      ! open(12, file = corr_file_name, status = "new")
      ! write(12, "(a, I7, a)") "Nz = ", Nz, ";   swil_corr_thin"
      ! close(12)

      open(12, file = swil_smeared_file_name, status = "new")
      write(12, "(a, I7, a)") "Nz = ", Nz, ";   swil_smeared"
      close(12)

      open(13, file = crosscorr_file_name, status = "new")
      write(13, "(a, I7, a)") "Nz = ", Nz, ";   swil_cross_corr"
      close(13)

      write(*,'(3(a,i0),2(a,f9.6))') ' su2_3D: Nx=',Nx,' Ny=',Ny,' Nz=',Nz,' beta_3D=',beta_3D,' hitsize=',hitsize
      write(*,'(6(a,i0),1(a,l1)  )') ' su2_3D: n_ther=',n_ther,' n_meas=',n_meas,' n_sepa=',n_sepa,' n_over=',n_over,' n_hit=',n_hit,' n_stout=',n_stout,' do_openmp=',do_openmp
      ! write(*,'(3(a,i0)          )') ' su2_3D: largest integer_sp is ',huge(1_sp),' and 2^31-1 is ',2147483647,' and largest box length is ',floor((dble(huge(1_sp))+1.0_dp)**0.33333)

#ifdef __GFORTRAN__
      write(*,'( (a)             )') ' su2_3D: this code has been compiled with GFORTRAN'
#elif  __INTEL_COMPILER
      write(*,'( (a)             )') ' su2_3D: this code has been compiled with IFORT'
#elif  __NAG__
      write(*,'( (a)             )') ' su2_3D: this code has been compiled with NAGFOR'
#elif __PGI__
      write(*,'( (a)             )') ' su2_3D: this code has been compiled with PGFORTRAN'
#endif

      call get_environment_variable('OMP_PLACES'   ,value=str_places)
      call get_environment_variable('OMP_PROC_BIND',value=str_procbind)
      call get_environment_variable('OMP_DYNAMIC'  ,value=str_dynamic)
      call get_environment_variable('KMP_AFFINITY' ,value=str_kmp)
      write(*,*) 'su2_3D: OMP_{PLACES,PROC_BIND,DYNAMIC}: '//trim(str_places)//' '//trim(str_procbind)//' '//trim(str_dynamic)
      write(*,*) 'su2_3D: KMP_AFFINITY: '//trim(str_kmp)

      call get_environment_variable('HOSTNAME',value=hostname)
      call date_and_time(date=hostdate,time=hosttime,zone=hostzone)
      write(*,*) 'su2_3D: host{name,date,time,zone}: '//trim(hostname)//' '//trim(hostdate)//' '//trim(hosttime)//' '//trim(hostzone)

      call system_clock(t_ini,clock_rate) !!! note: beg overall timing
      call su2_3D_init(U)

#ifdef _OPENMP
      n_procs=omp_get_num_procs()
      n_threads=omp_get_max_threads(); call omp_set_num_threads(n_threads)
#else
      n_procs=1
      n_threads=1
#endif
      write(*,'(3(a,i0)               )') ' su2_3D: henceforth this code runs with ',n_threads,' threads on ',n_procs,' virtual cores'

      do loop_ther=1,n_ther
         accrate_metro=0.0; accrate_ovlax=0.0
         call system_clock(t_new,clock_rate); t_min=huge(kind(t_min))
         do loop_sepa=1,n_sepa
            t_old=t_new
            call su2_3D_metropolis  (U,beta_3D,hitsize,n_hit,tmp_metro); accrate_metro=accrate_metro+tmp_metro/float(n_sepa)
            do loop_over=1,n_over
               call su2_3D_overrelax(U,beta_3D,              tmp_ovlax); accrate_ovlax=accrate_ovlax+tmp_ovlax/float(n_sepa*n_over)
            end do
            call system_clock(t_new,clock_rate); t_min=min(t_new-t_old,t_min)
         end do ! loop_sepa=1:n_sepa
         call su2_3D_backproject(U)
         hitsize=hitsize*sqrt(sqrt(accrate_metro/0.8))
         t_sec=real(dble(t_min)/dble(clock_rate),kind=sp)
         if (modulo(loop_ther,10)==0) then
            write(*,'(i8,a,2f7.4,a,f8.6,a,f8.6)') loop_ther-n_ther,' accrates=',accrate_metro,accrate_ovlax,' time=',t_sec,' hitsize= ',hitsize
         end if
      end do ! loop_ther=1:n_ther

      open(10, file = log_file_name, status = "old", position = "append")
      write(10, "(a, F7.4)") "hitsize after therm. = ", hitsize
      close(10)

      allocate(swil_smeared(1:size(smearlist)))
      allocate(swil_cross_corr(1:size(smearlist),1:size(smearlist),1:Nz))

      do loop_meas=1,n_meas
         accrate_metro=0.0; accrate_ovlax=0.0
         call system_clock(t_new,clock_rate); t_min=huge(kind(t_min))
         do loop_sepa=1,n_sepa
            t_old=t_new
            call su2_3D_metropolis  (U,beta_3D,hitsize,n_hit,tmp_metro); accrate_metro=accrate_metro+tmp_metro/float(n_sepa)
            do loop_over=1,n_over
               call su2_3D_overrelax(U,beta_3D              ,tmp_ovlax); accrate_ovlax=accrate_ovlax+tmp_ovlax/float(n_sepa*n_over)
            end do
            call system_clock(t_new,clock_rate); t_min=min(t_new-t_old,t_min)
         end do ! loop_sepa=1:n_sepa
         call su2_3D_backproject(U)
         ! call su2_2D_multistout(U,V,rho_2D,n_stout) !!! note: rho_2D<0.24, rho_3D<0.16, rho_4D<0.12
         call calc_swil_corr_mat_2D(U,smearlist,rho_2D,swil_smeared,swil_cross_corr)
         ! call calc_swil_corr_2D(U, swil_thin, swil_corr_thin)
         ! call calc_swil_corr_2D(V, swil_smth, swil_corr_smth)
         ! swil_thin=calc_swil(U)
         ! swil_smth=calc_swil(V)
         sopt_thin=calc_sopt(U)
         ! sopt_smth=calc_sopt(V)
         t_sec=real(dble(t_min)/dble(clock_rate),kind=sp)
         !!! order in timeseries_file: "s_wil    clover   smeared_s_wil  smeared_clover  accrate_metro  accrate_ovlax  time"
         open(10, file = timeseries_file_name, status = "old", position = "append")
         write(10,*) swil_smeared(1), sopt_thin, swil_smeared(size(smearlist)), sopt_smth, accrate_metro, accrate_ovlax, t_sec
         close(10) 
         ! open(11, file = corr_file_name, status = "old", position = "append")
         ! write(11,*) swil_corr_thin, swil_corr_smth
         ! close(11) 
         open(11, file = swil_smeared_file_name, status = "old", position = "append")
         write(11,*) swil_smeared
         close(11) 
         open(12, file = crosscorr_file_name, status = "old", position = "append")
         do z_loop = 1,Nz
            write(12,*) swil_cross_corr(:,:,z_loop)
         end do
         write(12,"(a)") " "
         close(12)
         if (modulo(loop_meas,10)==0) then
            write(*,'(i8,a,2f7.4,a,f8.6,a,2f9.6,a,2f9.6)') loop_meas,' accrates=',accrate_metro,accrate_ovlax,' time=',t_sec,' swil=',swil_thin,swil_smth,' sopt=',sopt_thin,sopt_smth
         end if
      end do ! loop_meas=1:n_meas

      write(*,'(3(a,i0),2(a,f9.6))') ' su2_3D: Nx=',Nx,' Ny=',Ny,' Nz=',Nz,' beta_3D=',beta_3D,' hitsize=',hitsize
      write(*,'(6(a,i0),1(a,l1)  )') ' su2_3D: n_ther=',n_ther,' n_meas=',n_meas,' n_sepa=',n_sepa,' n_over=',n_over,' n_hit=',n_hit,' n_stout=',n_stout,' do_openmp=',do_openmp
      call date_and_time(date=hostdate,time=hosttime,zone=hostzone)
      write(*,*) 'su2_3D: host{name,date,time,zone}: '//trim(hostname)//' '//trim(hostdate)//' '//trim(hosttime)//' '//trim(hostzone)

      call system_clock(t_fin,clock_rate) !!! note: end overall timing
      t_sec=real(dble(t_fin-t_ini)/dble(clock_rate),kind=sp)
      write(*,'(a,f10.2,f8.2,f6.2)') ' overall time in seconds/minutes/hours:',t_sec,t_sec/60.0,t_sec/3600.0
      write(*,*)
      open(10, file = log_file_name, status = "old", position = "append")
      write(10, "(a,f10.2,a,f8.2,a,f6.2)") "Overall time in seconds/minutes/hours: ",t_sec," / ",t_sec/60.0," / ",t_sec/3600.0

   end subroutine su2_3D_main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine su2_3D_init(U)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(out) :: U
      integer :: dir,x,y,z
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) SHARED(U)
      do z=1,Nz
      do y=1,Ny
      do x=1,Nx
      do dir=1,3
         U(:,dir,x,y,z)=cmplx([1.0,0.0],kind=sp)
      end do ! dir=1:3
      end do ! x=1:Nx
      end do ! y=1:Ny
      end do ! z=1:Nz
      !$OMP END PARALLEL DO
      write(*,*) 'su2_3D: initialization of U is cold'
   end subroutine su2_3D_init

#ifdef never
   subroutine myrandom_number_sp(rand_sp)
      implicit none
      real(kind=sp),dimension(:,:,:),intent(out) :: rand_sp
#ifdef __GFORTRAN__
      integer :: y,z
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) SHARED(rand_sp)
      do z=lbound(rand_sp,dim=3),ubound(rand_sp,dim=3)
      do y=lbound(rand_sp,dim=2),ubound(rand_sp,dim=2)
         call random_number(rand_sp(:,y,z))
      end do
      end do
      !$OMP END PARALLEL DO
#else
      integer(kind=dp),dimension(Ny,Nz),save :: seed_dp=0_dp
      integer(kind=dp) :: s_dp
      integer :: x,y,z,io_unit,io_stat
      real(kind=sp) :: r_sp
      if (all(seed_dp==0)) then
         open(newunit=io_unit,file='/dev/urandom',access='stream',form='unformatted',action='read',status='old',iostat=io_stat)
         if (io_stat==0) then !!! if OS provides a random number generator use it
            read(io_unit) seed_dp(:)
            seed_dp(:)=abs(seed_dp(:)) !!! note: avoid negative values as deilvered from "/dev/urandom"
            write(*,*) 'su2_3D: initialization of shared seed_dp through /dev/urandom gave',minval(seed_dp),maxval(seed_dp)
         else !!! otherwise use serial intrinsic function random_number to generate seed for each column
            do z=1,Nz
            do y=1,Ny
               call random_number(r_sp)
               seed_dp(y,z)=int(r_sp*huge(seed_dp),kind=dp)
            end do
            end do
            write(*,*) 'su2_3D: initialization of shared seed_dp through fortran-rand gave',minval(seed_dp),maxval(seed_dp)
         end if
         close(io_unit)
      end if
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(s_dp) SHARED(rand_sp,seed_dp)
      do z=lbound(rand_sp,dim=3),ubound(rand_sp,dim=3)
      do y=lbound(rand_sp,dim=2),ubound(rand_sp,dim=2)
         s_dp=seed_dp(y,z)
         do x=lbound(rand_sp,dim=1),ubound(rand_sp,dim=1)
            select case(2)
               case(1)
               s_dp=mod(1664525_dp*s_dp+1013904223,2147483647_dp)
               !rand_sp(x,y,z)=    float(s_dp)/2147483647.0     !!! note: 2^31-1=2147483647
               rand_sp(x,y,z)=sngl(dble(s_dp)/2147483647.0_dp) !!! note: 2^31-1=2147483647
               case(2)
               s_dp=iand(1103515245_dp*s_dp+12345,2147483647_dp) !!! note: iand(.,m-1) means mod(.,m) for m=2^integer
               !rand_sp(x,y,z)=    float(s_dp)/2147483648.0
               rand_sp(x,y,z)=sngl(dble(s_dp)/2147483648.0_dp)
            end select
         end do
         seed_dp(y,z)=s_dp
      end do
      end do
      !$OMP END PARALLEL DO
#endif
   end subroutine myrandom_number_sp
#endif

   pure function su2_full(a)
      !!! note: a must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a
      complex(kind=sp),dimension(2,2) :: su2_full
      su2_full(:,1)=a
      su2_full(:,2)=[-conjg(a(2)),conjg(a(1))]
   end function su2_full

   pure function su2_dagger(a)
      !!! note: a must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a
      complex(kind=sp),dimension(2) :: su2_dagger
      su2_dagger=[conjg(a(1)),-a(2)]
   end function su2_dagger

   pure function su2_normsqu(a)
      !!! note: a must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a
      real(kind=sp) :: su2_normsqu
      su2_normsqu=sum(realpart(a)**2+aimag(a)**2)
   end function su2_normsqu

   pure function su2_project(a)
      !!! note: a must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a
      complex(kind=sp),dimension(2) :: su2_project
      real(kind=sp) :: normsqu
      normsqu=su2_normsqu(a)
      su2_project=a/sqrt(normsqu) !!! note: specific to SU(2)
   end function su2_project

   pure function su2_realtrace(a)
      !!! note: a must be 2x1 column vector, and res is real scalar
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a
      real(kind=sp) :: su2_realtrace
      su2_realtrace=2.0*realpart(a(1))
   end function su2_realtrace

   pure function su2_mult_oo(a,b)
      !!! note: a,b,c must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a,b
      complex(kind=sp),dimension(2) :: su2_mult_oo
      su2_mult_oo=[a(1)*b(1)-conjg(a(2))*b(2),a(2)*b(1)+conjg(a(1))*b(2)]
   end function su2_mult_oo

   pure function su2_mult_do(a,b)
      !!! note: a,b,c must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a,b
      complex(kind=sp),dimension(2) :: su2_mult_do,a_dag
      a_dag=[conjg(a(1)),-a(2)]
      su2_mult_do=[a_dag(1)*b(1)-conjg(a_dag(2))*b(2),a_dag(2)*b(1)+conjg(a_dag(1))*b(2)]
   end function su2_mult_do

   pure function su2_mult_od(a,b)
      !!! note: a,b,c must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a,b
      complex(kind=sp),dimension(2) :: su2_mult_od,b_dag
      b_dag=[conjg(b(1)),-b(2)]
      su2_mult_od=[a(1)*b_dag(1)-conjg(a(2))*b_dag(2),a(2)*b_dag(1)+conjg(a(1))*b_dag(2)]
   end function su2_mult_od

   pure function su2_mult_ooo(a,b,c)
      !!! note: a,b,c must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a,b,c
      complex(kind=sp),dimension(2) :: su2_mult_ooo,tmp
      tmp=su2_mult_oo(b,c)
      su2_mult_ooo=su2_mult_oo(a,tmp)
   end function su2_mult_ooo

   pure function su2_mult_doo(a,b,c)
      !!! note: a,b,c must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a,b,c
      complex(kind=sp),dimension(2) :: su2_mult_doo,a_dag,tmp
      a_dag=[conjg(a(1)),-a(2)] !!! note: su2_dagger(a)
      tmp=su2_mult_oo(b,c)
      su2_mult_doo=su2_mult_oo(a_dag,tmp)
   end function su2_mult_doo

   pure function su2_mult_odo(a,b,c)
      !!! note: a,b,c must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a,b,c
      complex(kind=sp),dimension(2) :: su2_mult_odo,b_dag,tmp
      b_dag=[conjg(b(1)),-b(2)] !!! note: su2_dagger(b)
      tmp=su2_mult_oo(b_dag,c)
      su2_mult_odo=su2_mult_oo(a,tmp)
   end function su2_mult_odo

   pure function su2_mult_ood(a,b,c)
      !!! note: a,b,c must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a,b,c
      complex(kind=sp),dimension(2) :: su2_mult_ood,c_dag,tmp
      c_dag=[conjg(c(1)),-c(2)] !!! note: su2_dagger(c)
      tmp=su2_mult_oo(b,c_dag)
      su2_mult_ood=su2_mult_oo(a,tmp)
   end function su2_mult_ood

   pure function su2_mult_oodd(a,b,c,d)
      !!! note: a,b,c must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a,b,c,d
      complex(kind=sp),dimension(2) :: su2_mult_oodd,c_dag,d_dag,tmp,tmq
      c_dag=[conjg(c(1)),-c(2)] !!! note: su2_dagger(c)
      d_dag=[conjg(d(1)),-d(2)] !!! note: su2_dagger(d)
      tmp=su2_mult_oo(c_dag,d_dag)
      tmq=su2_mult_oo(b,tmp)
      su2_mult_oodd=su2_mult_oo(a,tmq)
   end function su2_mult_oodd

   pure function su2_mult_oddo(a,b,c,d)
      !!! note: a,b,c must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a,b,c,d
      complex(kind=sp),dimension(2) :: su2_mult_oddo,b_dag,c_dag,tmp,tmq
      b_dag=[conjg(b(1)),-b(2)] !!! note: su2_dagger(b)
      c_dag=[conjg(c(1)),-c(2)] !!! note: su2_dagger(c)
      tmp=su2_mult_oo(c_dag,d)
      tmq=su2_mult_oo(b_dag,tmp)
      su2_mult_oddo=su2_mult_oo(a,tmq)
   end function su2_mult_oddo

   pure function su2_mult_ddoo(a,b,c,d)
      !!! note: a,b,c must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a,b,c,d
      complex(kind=sp),dimension(2) :: su2_mult_ddoo,a_dag,b_dag,tmp,tmq
      a_dag=[conjg(a(1)),-a(2)] !!! note: su2_dagger(a)
      b_dag=[conjg(b(1)),-b(2)] !!! note: su2_dagger(b)
      tmp=su2_mult_oo(c,d)
      tmq=su2_mult_oo(b_dag,tmp)
      su2_mult_ddoo=su2_mult_oo(a_dag,tmq)
   end function su2_mult_ddoo

   pure function su2_mult_dood(a,b,c,d)
      !!! note: a,b,c must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a,b,c,d
      complex(kind=sp),dimension(2) :: su2_mult_dood,a_dag,d_dag,tmp,tmq
      a_dag=[conjg(a(1)),-a(2)] !!! note: su2_dagger(a)
      d_dag=[conjg(d(1)),-d(2)] !!! note: su2_dagger(d)
      tmp=su2_mult_oo(c,d_dag)
      tmq=su2_mult_oo(b,tmp)
      su2_mult_dood=su2_mult_oo(a_dag,tmq)
   end function su2_mult_dood

   pure function su2_random(hitsize,r)
      implicit none
      real(kind=sp),intent(in) :: hitsize,r(3)
      real(kind=sp) :: s(3),fac_hitsize
      complex(kind=sp),dimension(2) :: su2_random
      s=r-0.5; s=(hitsize/sqrt(sum(s**2)))*s
      fac_hitsize=sqrt(1.0-hitsize**2)
      !!! res=[complex(fac_hitsize,+s(3)),complex(+s(2),s(1)); ...
      !!!      complex(-s(2),s(1)),complex(fac_hitsize,-s(3))];
      su2_random=[cmplx(fac_hitsize,s(3),kind=sp),cmplx(-s(2),s(1),kind=sp)]
   end function su2_random

   pure function su2_inverse(a)
      !!! note: a must be 2x1 column vector, and res is also (2x1)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a
      complex(kind=sp),dimension(2) :: su2_inverse
      !!! note: works also if argument A is MULTIPLE of SU(2) matrix
      su2_inverse=[conjg(a(1)),-a(2)]/su2_normsqu(a)
   end function su2_inverse

   pure function su2_expm(a)
      implicit none
      complex(kind=sp),dimension(2),intent(in) :: a
      complex(kind=sp),dimension(2) :: su2_expm,num,den,lastterm
      real(kind=sp) :: fac,lastnorm
      integer :: n
      !!! note: argument A is anti-hermitean and traceless, NOT in SU(2)
      lastterm=[cmplx(1.0,kind=sp),cmplx(0.0,kind=sp)]
      num=lastterm !!! note: to become expm(+A/8), take num^4 in the end
      den=lastterm !!! note: to become expm(-A/8), take den^4 in the end
      fac=1.0
      do n=1,100
         lastterm=su2_mult_oo(lastterm,a)/float(8*n)
         lastnorm=sqrt(su2_normsqu(lastterm))
         fac=-fac
         num=num+1.0*lastterm
         den=den+fac*lastterm
         if (lastnorm.lt.1e-18) exit
      end do
      den=su2_inverse(den) !!! note: den on RHS is multiple of SU(2)
      su2_expm=su2_mult_oo(num,den)
      su2_expm=su2_mult_oo(su2_expm,su2_expm)
      su2_expm=su2_mult_oo(su2_expm,su2_expm)
      su2_expm=su2_project(su2_expm)
   end function su2_expm

   pure function su2_3D_staple(U,dir,x,y,z)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(in) :: U
      integer,intent(in) :: dir,x,y,z
      complex(kind=sp),dimension(2) :: su2_3D_staple
      integer :: x_min,x_plu,y_min,y_plu,z_min,z_plu
      x_min=modulo(x-2,Nx)+1; x_plu=modulo(x,Nx)+1;
      y_min=modulo(y-2,Ny)+1; y_plu=modulo(y,Ny)+1;
      z_min=modulo(z-2,Nz)+1; z_plu=modulo(z,Nz)+1;
#define oldway
#ifdef oldway
      select case(dir)
         case(1)
         su2_3D_staple=su2_mult_ood(U(:,2,x,y    ,z),U(:,1,x,y_plu,z),U(:,2,x_plu,y    ,z)) &
                        +su2_mult_doo(U(:,2,x,y_min,z),U(:,1,x,y_min,z),U(:,2,x_plu,y_min,z)) &
                        +su2_mult_ood(U(:,3,x,y,z    ),U(:,1,x,y,z_plu),U(:,3,x_plu,y,z    )) &
                        +su2_mult_doo(U(:,3,x,y,z_min),U(:,1,x,y,z_min),U(:,3,x_plu,y,z_min))
         case(2)
         su2_3D_staple=su2_mult_ood(U(:,3,x,y,z    ),U(:,2,x,y,z_plu),U(:,3,x,y_plu,z    )) &
                        +su2_mult_doo(U(:,3,x,y,z_min),U(:,2,x,y,z_min),U(:,3,x,y_plu,z_min)) &
                        +su2_mult_ood(U(:,1,x    ,y,z),U(:,2,x_plu,y,z),U(:,1,x    ,y_plu,z)) &
                        +su2_mult_doo(U(:,1,x_min,y,z),U(:,2,x_min,y,z),U(:,1,x_min,y_plu,z))
         case(3)
         su2_3D_staple=su2_mult_ood(U(:,1,x    ,y,z),U(:,3,x_plu,y,z),U(:,1,x    ,y,z_plu)) &
                        +su2_mult_doo(U(:,1,x_min,y,z),U(:,3,x_min,y,z),U(:,1,x_min,y,z_plu)) &
                        +su2_mult_ood(U(:,2,x,y    ,z),U(:,3,x,y_plu,z),U(:,2,x,y    ,z_plu)) &
                        +su2_mult_doo(U(:,2,x,y_min,z),U(:,3,x,y_min,z),U(:,2,x,y_min,z_plu))
      end select
#else
      select case(dir)
         case(1)
         su2_3D_staple=su2_project(su2_mult_ood(U(:,2,x,y    ,z),U(:,1,x,y_plu,z),U(:,2,x_plu,y    ,z))) &
                        +su2_project(su2_mult_doo(U(:,2,x,y_min,z),U(:,1,x,y_min,z),U(:,2,x_plu,y_min,z))) &
                        +su2_project(su2_mult_ood(U(:,3,x,y,z    ),U(:,1,x,y,z_plu),U(:,3,x_plu,y,z    ))) &
                        +su2_project(su2_mult_doo(U(:,3,x,y,z_min),U(:,1,x,y,z_min),U(:,3,x_plu,y,z_min)))
         case(2)
         su2_3D_staple=su2_project(su2_mult_ood(U(:,3,x,y,z    ),U(:,2,x,y,z_plu),U(:,3,x,y_plu,z    ))) &
                        +su2_project(su2_mult_doo(U(:,3,x,y,z_min),U(:,2,x,y,z_min),U(:,3,x,y_plu,z_min))) &
                        +su2_project(su2_mult_ood(U(:,1,x    ,y,z),U(:,2,x_plu,y,z),U(:,1,x    ,y_plu,z))) &
                        +su2_project(su2_mult_doo(U(:,1,x_min,y,z),U(:,2,x_min,y,z),U(:,1,x_min,y_plu,z)))
         case(3)
         su2_3D_staple=su2_project(su2_mult_ood(U(:,1,x    ,y,z),U(:,3,x_plu,y,z),U(:,1,x    ,y,z_plu))) &
                        +su2_project(su2_mult_doo(U(:,1,x_min,y,z),U(:,3,x_min,y,z),U(:,1,x_min,y,z_plu))) &
                        +su2_project(su2_mult_ood(U(:,2,x,y    ,z),U(:,3,x,y_plu,z),U(:,2,x,y    ,z_plu))) &
                        +su2_project(su2_mult_doo(U(:,2,x,y_min,z),U(:,3,x,y_min,z),U(:,2,x,y_min,z_plu)))
      end select
#endif
      !!! note: staple is not in SU(2) but [specifically for SU(2)] a multiple of an SU(2) element, and
      !!! this is why the one-column representation can still be used, except that it is not normalized
   end function su2_3D_staple

   pure function su2_2D_staple(U,dir,x,y,z)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(in) :: U
      integer,intent(in) :: dir,x,y,z
      complex(kind=sp),dimension(2) :: su2_2D_staple
      integer :: x_min,x_plu,y_min,y_plu
      x_min=modulo(x-2,Nx)+1; x_plu=modulo(x,Nx)+1;
      y_min=modulo(y-2,Ny)+1; y_plu=modulo(y,Ny)+1;
#ifdef oldway
      select case(dir)
         case(1)
         su2_2D_staple=su2_mult_ood(U(:,2,x,y    ,z),U(:,1,x,y_plu,z),U(:,2,x_plu,y    ,z)) &
                        +su2_mult_doo(U(:,2,x,y_min,z),U(:,1,x,y_min,z),U(:,2,x_plu,y_min,z))
         case(2)
         su2_2D_staple=su2_mult_ood(U(:,1,x    ,y,z),U(:,2,x_plu,y,z),U(:,1,x    ,y_plu,z)) &
                        +su2_mult_doo(U(:,1,x_min,y,z),U(:,2,x_min,y,z),U(:,1,x_min,y_plu,z))
      end select
#else
      select case(dir)
         case(1)
         su2_2D_staple=su2_project(su2_mult_ood(U(:,2,x,y    ,z),U(:,1,x,y_plu,z),U(:,2,x_plu,y    ,z))) &
                        +su2_project(su2_mult_doo(U(:,2,x,y_min,z),U(:,1,x,y_min,z),U(:,2,x_plu,y_min,z)))
         case(2)
         su2_2D_staple=su2_project(su2_mult_ood(U(:,1,x    ,y,z),U(:,2,x_plu,y,z),U(:,1,x    ,y_plu,z))) &
                        +su2_project(su2_mult_doo(U(:,1,x_min,y,z),U(:,2,x_min,y,z),U(:,1,x_min,y_plu,z)))
      end select
#endif
      !!! note: staple is not in SU(2) but [specifically for SU(2)] a multiple of an SU(2) element, and
      !!! this is why the one-column representation can still be used, except that it is not normalized
   end function su2_2D_staple

   subroutine su2_3D_multistout(U,V,rho_3D,n_stout)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent( in) :: U
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(out) :: V
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz) :: V_tmp
      real(kind=sp),intent(in) :: rho_3D
      integer,intent(in) :: n_stout
      integer :: n
      call su2_3D_confcopy(U,V)
      do n=1,n_stout
         call su2_3D_confcopy(V,V_tmp)
         call su2_3D_stout(V_tmp,V,rho_3D) !!! note: must be smaller than 1/6
      end do
   end subroutine su2_3D_multistout

   subroutine su2_2D_multistout(U,V,rho_2D,n_stout)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent( in) :: U
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(out) :: V
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz) :: V_tmp
      real(kind=sp),intent(in) :: rho_2D
      integer,intent(in) :: n_stout
      integer :: n
      call su2_3D_confcopy(U,V)
      do n=1,n_stout
         call su2_3D_confcopy(V,V_tmp)
         call su2_2D_stout(V_tmp,V,rho_2D) !!! note: must be smaller than 1/4
      end do
   end subroutine su2_2D_multistout

   subroutine su2_3D_stout(V_tmp,V,rho_3D)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent( in) :: V_tmp
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(out) :: V
      real(kind=sp),intent(in) :: rho_3D
      complex(kind=sp),dimension(2) :: tmp
      integer :: dir,x,y,z
      if (rho_3D>1.0/6.0) then; write(*,*) 'su2_3D_stout: rho_3D<1/6 is violated'; stop; end if
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(tmp) FIRSTPRIVATE(rho_3D) SHARED(V_tmp,V) SCHEDULE(static)
      do z=1,Nz
      do y=1,Ny
      do x=1,Nx
      do dir=1,3
         tmp=su2_3D_staple(V_tmp,dir,x,y,z)
         tmp=su2_mult_od(tmp,V_tmp(:,dir,x,y,z))
         tmp=0.5*(tmp-su2_dagger(tmp))           !!! note: make anti-herm
         tmp(1)=cmplx(0.0,aimag(tmp(1)),kind=sp) !!! note: make traceless
         tmp=su2_expm(rho_3D*tmp)
         tmp=su2_mult_oo(tmp,V_tmp(:,dir,x,y,z))
         V(:,dir,x,y,z)=tmp
      end do ! dir=1:3
      end do ! x=1:Nx
      end do ! y=1:Ny
      end do ! z=1:Nz
      !$OMP END PARALLEL DO
   end subroutine su2_3D_stout

   subroutine su2_2D_stout(V_tmp,V,rho_2D)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent( in) :: V_tmp
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(out) :: V
      real(kind=sp),intent(in) :: rho_2D
      complex(kind=sp),dimension(2) :: tmp
      integer :: dir,x,y,z
      if (rho_2D>1.0/4.0) then; write(*,*) 'su2_2D_stout: rho_2D<1/4 is violated'; stop; end if
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(tmp) FIRSTPRIVATE(rho_2D) SHARED(V_tmp,V) SCHEDULE(static)
      do z=1,Nz
      do y=1,Ny
      do x=1,Nx
         do dir=1,2
            tmp=su2_2D_staple(V_tmp,dir,x,y,z)
            tmp=su2_mult_od(tmp,V_tmp(:,dir,x,y,z))
            tmp=0.5*(tmp-su2_dagger(tmp))           !!! note: make anti-herm
            tmp(1)=cmplx(0.0,aimag(tmp(1)),kind=sp) !!! note: make traceless
            tmp=su2_expm(rho_2D*tmp)
            tmp=su2_mult_oo(tmp,V_tmp(:,dir,x,y,z))
            V(:,dir,x,y,z)=tmp
         end do ! dir=1:2
         ! V(:,3,x,y,z) = V_tmp(:,3,x,y,z)  !!! no smearing in dir=3
      end do ! x=1:Nx
      end do ! y=1:Ny
      end do ! z=1:Nz
      !$OMP END PARALLEL DO
   end subroutine su2_2D_stout

   subroutine su2_3D_confcopy(V,W)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent( in) :: V
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(out) :: W
      integer :: x,y,z
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) SHARED(V,W) SCHEDULE(static)
      do z=1,Nz
      do y=1,Ny
      do x=1,Nx
         W(:,:,x,y,z)=V(:,:,x,y,z)
      end do ! x=1:Nx
      end do ! y=1:Ny
      end do ! z=1:Nz
      !$OMP END PARALLEL DO
   end subroutine su2_3D_confcopy

   subroutine su2_3D_metropolis(U,beta_3D,hitsize,n_hit,accrate)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(inout) :: U
      real(kind=sp)   ,intent( in) :: beta_3D,hitsize
      real(kind=sp)   ,intent(out) :: accrate
      integer         ,intent( in) :: n_hit
      complex(kind=sp),dimension(2) :: oldlink,newlink,staple_dag
      real(kind=sp),dimension(4,n_hit) :: r_sp
      real(kind=sp) :: s_old,s_new,boltz
      integer :: x,y,z,trip,dir,count_metro,loop_hit
      if (hitsize>1.0) then
         write(*,*) 'su2_3D_metropolis: hitsize<1.0 is violated by hitsize=',hitsize; stop
      end if
      count_metro=0
      do dir=1,3
      do trip=1,2
         !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(oldlink,newlink,staple_dag,s_old,s_new,boltz,r_sp) FIRSTPRIVATE(beta_3D,hitsize,n_hit,dir,trip) SHARED(U) REDUCTION(+:count_metro) SCHEDULE(static)
         do z=1,Nz
         do y=1,Ny
         do x=1+modulo(y+z+trip,2),Nx,2
            oldlink(:)=U(:,dir,x,y,z)
            staple_dag=su2_dagger(su2_3D_staple(U,dir,x,y,z)) !!! note: is (2x1) and NOT in SU(2)
            select case(2)
               case(1);    s_old=1.0-su2_realtrace(su2_mult_oo(oldlink,staple_dag))/2.0
               case(2);    s_old=1.0-realpart(oldlink(1)*staple_dag(1)-conjg(oldlink(2))*staple_dag(2))
            end select
            call random_number(r_sp)
            do loop_hit=1,n_hit
               newlink=su2_mult_oo(su2_random(hitsize,r_sp(1:3,loop_hit)),oldlink)
               newlink=su2_project(newlink)
               select case(2)
                  case(1); s_new=1.0-su2_realtrace(su2_mult_oo(newlink,staple_dag))/2.0
                  case(2); s_new=1.0-realpart(newlink(1)*staple_dag(1)-conjg(newlink(2))*staple_dag(2))
               end select
               boltz=exp(-beta_3D*(s_new-s_old))
               if (r_sp(4,loop_hit)<boltz) then
                  oldlink=newlink
                  s_old=s_new
                  count_metro=count_metro+1
               end if
            end do ! loop_hit=1:n_hit
            U(:,dir,x,y,z)=oldlink
         end do ! x=1:Nx
         end do ! y=1:Ny
         end do ! z=1:Nz
         !$OMP END PARALLEL DO
      end do ! trip=1:2
      end do ! dir=1:3
      accrate=real(dble(count_metro)/dble(n_hit*3*Nx*Ny*Nz),kind=sp)
   end subroutine su2_3D_metropolis

   subroutine su2_3D_overrelax(U,beta_3D,accrate)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(inout) :: U
      real(kind=sp),intent( in) :: beta_3D
      real(kind=sp),intent(out) :: accrate
      complex(kind=sp),dimension(2) :: oldlink,newlink,staple,staple_proj
      real(kind=sp) :: r_sp,boltz
      integer :: x,y,z,dir,trip,count_ovlax
      count_ovlax=0
      do dir=1,3
      do trip=1,2
         !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(oldlink,newlink,staple,staple_proj,boltz,r_sp) FIRSTPRIVATE(beta_3D,dir,trip) SHARED(U) REDUCTION(+:count_ovlax) SCHEDULE(static)
         do z=1,Nz
         do y=1,Ny
         do x=1+modulo(y+z+trip,2),Nx,2
            staple=su2_3D_staple(U,dir,x,y,z) !!! note: is (2x1) and not in SU(2)
            staple_proj=su2_project(staple)
            oldlink=U(:,dir,x,y,z)
            newlink=su2_mult_odo(staple_proj,oldlink,staple_proj)
            newlink=su2_project(newlink)
            boltz=exp(beta_3D*su2_realtrace(su2_mult_od(newlink-oldlink,staple))/2.0)
            call random_number(r_sp)
            if (r_sp<boltz) then
               U(:,dir,x,y,z)=newlink
               count_ovlax=count_ovlax+1
            end if
         end do ! x=1:Nx
         end do ! y=1:Ny
         end do ! z=1:Nz
         !$OMP END PARALLEL DO
      end do ! trip=1:2
      end do ! dir=1:3
      accrate=real(dble(count_ovlax)/dble(3*Nx*Ny*Nz),kind=sp)
   end subroutine su2_3D_overrelax

   subroutine su2_3D_backproject(U)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(inout) :: U
      integer :: x,y,z,dir
      do dir=1,3
         !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) FIRSTPRIVATE(dir) SHARED(U) SCHEDULE(static)
         do z=1,Nz
         do y=1,Ny
         do x=1,Nx
            U(:,dir,x,y,z)=su2_project(U(:,dir,x,y,z))
         end do ! x=1:Nx
         end do ! y=1:Ny
         end do ! z=1:Nz
         !$OMP END PARALLEL DO
      end do ! dir=1:3
   end subroutine su2_3D_backproject

   function calc_swil(U)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(in) :: U
      real(kind=sp) :: calc_swil
      real(kind=dp) :: s_dp
      integer :: x,y,z,x_plu,y_plu,z_plu
      s_dp=0.0_dp
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(x_plu,y_plu,z_plu) SHARED(U) REDUCTION(+:s_dp) SCHEDULE(static)
      do z=1,Nz
      do y=1,Ny; y_plu=modulo(y,Ny)+1; z_plu=modulo(z,Nz)+1
      do x=1,Nx; x_plu=modulo(x,Nx)+1
         s_dp=s_dp-su2_realtrace(su2_mult_oodd(U(:,1,x,y,z),U(:,2,x_plu,y,z),U(:,1,x,y_plu,z),U(:,2,x,y,z))) &
                  -su2_realtrace(su2_mult_oodd(U(:,1,x,y,z),U(:,3,x_plu,y,z),U(:,1,x,y,z_plu),U(:,3,x,y,z))) &
                  -su2_realtrace(su2_mult_oodd(U(:,2,x,y,z),U(:,3,x,y_plu,z),U(:,2,x,y,z_plu),U(:,3,x,y,z)))
      end do ! x=1:Nx
      end do ! y=1:Ny
      end do ! z=1:Nz
      !$OMP END PARALLEL DO
      calc_swil=real(dble(1)+s_dp/dble(2*3*Nx*Ny*Nz),kind=sp)
   end function calc_swil

   subroutine calc_swil_corr_2D(U,s_dp,corr)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(in) :: U
      real(kind=dp), intent(out) :: s_dp, corr(1:Nz)
      real(kind=dp) :: retrp, s_slice(Nz) 
      integer :: x, y, z, x_plu, y_plu, z_plu
      !$OMP PARALLEL DO DEFAULT(none) PRIVATE(x_plu,y_plu,z_plu,retrp) SHARED(U, s_slice) SCHEDULE(static)
      do z=1,Nz; z_plu=modulo(z,Nz)+1
         retrp = 0.0
         do y=1,Ny; y_plu=modulo(y,Ny)+1 
         do x=1,Nx; x_plu=modulo(x,Nx)+1 
            retrp=retrp+su2_realtrace(su2_mult_oodd(U(:,1,x,y,z),U(:,2,x_plu,y,z),U(:,1,x,y_plu,z),U(:,2,x,y,z)))
         end do ! x=1:Nx
         end do ! y=1:Ny
         s_slice(z) = 1-retrp/(2*Nx*Ny)
      end do ! z=1:Nz
      !$OMP END PARALLEL DO
      do z=1,Nz
         corr(z) = sum(s_slice * cshift(s_slice,z))/Nz
      end do
      s_dp = sum(s_slice)/dble(Nz)
   end subroutine calc_swil_corr_2D

   subroutine calc_swil_corr_mat_2D(U,smearlist,rho_2D,s_smeared_dp,corr)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(in) :: U
      integer, intent(in) :: smearlist(:)
      real(kind=sp), intent(in)  :: rho_2D
      real(kind=dp), intent(out) :: s_smeared_dp(1:size(smearlist)), corr(1:size(smearlist),1:size(smearlist),1:Nz)
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz) :: V, V_tmp
      real(kind=dp) :: retrp, s_slice(1:Nz,1:size(smearlist))
      integer :: x, y, z, x_plu, y_plu, z_plu, smearhelp(1:size(smearlist)), i, smear_loop, op1, op2
      retrp = 0.0_dp
      smearhelp(1) = 0
      do i = 2,size(smearlist)
         smearhelp(i) = smearlist(i) - smearlist(i-1)
      end do
      call su2_3D_confcopy(U,V)
      !!! It's faster to not parallelize the smear_loop, only again the z-loop inside (for up to 5 smearing, up to Nz = 32)
      do smear_loop=1,size(smearhelp); call su2_2D_multistout(V,V_tmp,rho_2D,smearhelp(smear_loop))
         !$OMP PARALLEL DO DEFAULT(none) PRIVATE(x_plu,y_plu,z_plu,retrp) SHARED(V_tmp,s_slice,smear_loop) SCHEDULE(static)
         do z=1,Nz; z_plu=modulo(z,Nz)+1
            retrp = 0.0_dp
            do y=1,Ny; y_plu=modulo(y,Ny)+1
            do x=1,Nx; x_plu=modulo(x,Nx)+1
               retrp=retrp+su2_realtrace(su2_mult_oodd(V_tmp(:,1,x,y,z),V_tmp(:,2,x_plu,y,z),V_tmp(:,1,x,y_plu,z),V_tmp(:,2,x,y,z)))
            end do ! x=1:Nx
            end do ! y=1:Ny
            s_slice(z,smear_loop) = 1-retrp/(2*Nx*Ny)
         end do ! z=1:Nz
         !$OMP END PARALLEL DO
         call su2_3D_confcopy(V_tmp,V)
      end do ! smear_loop=1:size(smearhelp)
      do op1 = 1,size(smearlist)
      do op2 = 1,size(smearlist)
      do z=1,Nz
         corr(op1,op2,z) = sum(s_slice(:,op1) * cshift(s_slice(:,op2),z))/Nz
      end do ! z=0:Nz-1
      end do ! op2=1:size(smearlist)
      end do ! op1=1:size(smearlist)
      do i = 1,size(smearlist)
         s_smeared_dp(i) = sum(s_slice(:,i))/dble(Nz)
      end do
   end subroutine calc_swil_corr_mat_2D

   ! subroutine calc_swil_corr_mat_2D(U,smearlist,rho_2D,s_dp,s_smeared_dp,corr)
   !    implicit none
   !    complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(in) :: U
   !    integer, intent(in) :: smearlist(:)
   !    real(kind=sp), intent(in)  :: rho_2D
   !    real(kind=dp), intent(out) :: s_dp, s_smeared_dp, corr(1:size(smearlist),1:size(smearlist),1:Nz)
   !    complex(kind=sp),dimension(2,3,Nx,Ny,Nz) :: V, V_tmp
   !    real(kind=dp) :: retrp, s_slice(1:Nz,1:size(smearlist))
   !    integer :: x, y, z, x_plu, y_plu, z_plu, smearhelp(1:size(smearlist)), i, smear_loop, op1, op2
   !    retrp = 0.0_dp
   !    smearhelp(1) = 0
   !    do i = 2,size(smearlist)
   !       smearhelp(i) = smearlist(i) - smearlist(i-1)
   !    end do
   !    call su2_3D_confcopy(U,V)
   !    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(x_plu,y_plu,z_plu,retrp,V_tmp) SHARED(s_slice,smearlist,rho_2D,V) SCHEDULE(static)
   !    do smear_loop=1,size(smearlist); call su2_2D_multistout(V,V_tmp,rho_2D,smearlist(smear_loop))
   !       do z=1,Nz; z_plu=modulo(z,Nz)+1
   !          retrp = 0.0_dp
   !          do y=1,Ny; y_plu=modulo(y,Ny)+1
   !          do x=1,Nx; x_plu=modulo(x,Nx)+1
   !             retrp=retrp+su2_realtrace(su2_mult_oodd(V_tmp(:,1,x,y,z),V_tmp(:,2,x_plu,y,z),V_tmp(:,1,x,y_plu,z),V_tmp(:,2,x,y,z)))
   !          end do ! x=1:Nx
   !          end do ! y=1:Ny
   !          s_slice(z,smear_loop) = 1-retrp/(2*Nx*Ny)
   !       end do ! z=1:Nz
   !       ! call su2_3D_confcopy(V_tmp,V)
   !    end do ! smear_loop=1:size(smearhelp)
   !    !$OMP END PARALLEL DO
   !    do op1 = 1,size(smearlist)
   !    do op2 = 1,size(smearlist)
   !    do z=1,Nz
   !       corr(op1,op2,z) = sum(s_slice(:,op1) * cshift(s_slice(:,op2),z))/Nz
   !    end do ! z=0:Nz-1
   !    end do ! op2=1:size(smearlist)
   !    end do ! op1=1:size(smearlist)
   !    s_dp         = sum(s_slice(:,1))/dble(Nz)
   !    s_smeared_dp = sum(s_slice(:,size(smearlist)))/dble(Nz)
   ! end subroutine calc_swil_corr_mat_2D

   function calc_sopt(U)
      implicit none
      complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(in) :: U
      real(kind=sp) :: calc_sopt
      complex(kind=sp) :: tmq(2),tmr(2,2)
      real(kind=dp) :: s_dp
      integer :: x,y,z,x_min,y_min,z_min,x_plu,y_plu,z_plu
      s_dp=0.0_dp
      !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(x_min,x_plu,y_min,y_plu,z_min,z_plu,tmq,tmr) SHARED(U) REDUCTION(+:s_dp) SCHEDULE(static)
      do z=1,Nz
      do y=1,Ny; y_min=modulo(y-2,Ny)+1; y_plu=modulo(y,Ny)+1; z_min=modulo(z-2,Nz)+1; z_plu=modulo(z,Nz)+1
      do x=1,Nx; x_min=modulo(x-2,Nx)+1; x_plu=modulo(x,Nx)+1
         tmq=su2_mult_oodd(U(:,1,x,y,z),U(:,2,x_plu,y,z),U(:,1,x,y_plu,z),U(:,2,x,y,z)) &
            +su2_mult_oddo(U(:,2,x,y,z),U(:,1,x_min,y_plu,z),U(:,2,x_min,y,z),U(:,1,x_min,y,z)) &
            +su2_mult_ddoo(U(:,1,x_min,y,z),U(:,2,x_min,y_min,z),U(:,1,x_min,y_min,z),U(:,2,x,y_min,z)) &
            +su2_mult_dood(U(:,2,x,y_min,z),U(:,1,x,y_min,z),U(:,2,x_plu,y_min,z),U(:,1,x,y,z))
         tmr=su2_full(tmq)/cmplx(0.0,4.0,kind=sp); tmr=0.5*(tmr+conjg(transpose(tmr))); tmr(1,1)=0.5*(tmr(1,1)-tmr(2,2)); tmr(2,2)=-tmr(1,1) !!! note: make traceless
         s_dp=s_dp+realpart(sum(tmr(1,:)*tmr(:,1))+sum(tmr(2,:)*tmr(:,2)))/(2.0*2.0)
         tmq=su2_mult_oodd(U(:,1,x,y,z),U(:,3,x_plu,y,z),U(:,1,x,y,z_plu),U(:,3,x,y,z)) &
            +su2_mult_oddo(U(:,3,x,y,z),U(:,1,x_min,y,z_plu),U(:,3,x_min,y,z),U(:,1,x_min,y,z)) &
            +su2_mult_ddoo(U(:,1,x_min,y,z),U(:,3,x_min,y,z_min),U(:,1,x_min,y,z_min),U(:,3,x,y,z_min)) &
            +su2_mult_dood(U(:,3,x,y,z_min),U(:,1,x,y,z_min),U(:,3,x_plu,y,z_min),U(:,1,x,y,z))
         tmr=su2_full(tmq)/cmplx(0.0,4.0,kind=sp); tmr=0.5*(tmr+conjg(transpose(tmr))); tmr(1,1)=0.5*(tmr(1,1)-tmr(2,2)); tmr(2,2)=-tmr(1,1) !!! note: make traceless
         s_dp=s_dp+realpart(sum(tmr(1,:)*tmr(:,1))+sum(tmr(2,:)*tmr(:,2)))/(2.0*2.0)
         tmq=su2_mult_oodd(U(:,2,x,y,z),U(:,3,x,y_plu,z),U(:,2,x,y,z_plu),U(:,3,x,y,z)) &
            +su2_mult_oddo(U(:,3,x,y,z),U(:,2,x,y_min,z_plu),U(:,3,x,y_min,z),U(:,2,x,y_min,z)) &
            +su2_mult_ddoo(U(:,2,x,y_min,z),U(:,3,x,y_min,z_min),U(:,2,x,y_min,z_min),U(:,3,x,y,z_min)) &
            +su2_mult_dood(U(:,3,x,y,z_min),U(:,2,x,y,z_min),U(:,3,x,y_plu,z_min),U(:,2,x,y,z))
         tmr=su2_full(tmq)/cmplx(0.0,4.0,kind=sp); tmr=0.5*(tmr+conjg(transpose(tmr))); tmr(1,1)=0.5*(tmr(1,1)-tmr(2,2)); tmr(2,2)=-tmr(1,1) !!! note: make traceless
         s_dp=s_dp+realpart(sum(tmr(1,:)*tmr(:,1))+sum(tmr(2,:)*tmr(:,2)))/(2.0*2.0)
      end do ! x=1:Nx
      end do ! y=1:Ny
      end do ! z=1:Nz
      !$OMP END PARALLEL DO
      calc_sopt=real(s_dp/dble(3*Nx*Ny*Nz),kind=sp)
   end function calc_sopt

   ! subroutine calc_swil_sopt_corr_2D(U,s_dp,corr)
   !    implicit none
   !    complex(kind=sp),dimension(2,3,Nx,Ny,Nz),intent(in) :: U
   !    real(kind=dp), intent(out) :: s_dp, corr(Nz)
   !    real(kind=dp) :: s_slice(Nz), s_dp_old
   !    integer :: x, y, z, x_plu, y_plu, z_plu
   !    s_dp = 0.0_dp
   !    s_dp_old = 0.0_dp
   !    !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(none) PRIVATE(x_plu,y_plu,z_plu) SHARED(U, s_slice) REDUCTION(+:s_dp, s_dp_old) SCHEDULE(static)
   !    do z=1,Nz; 
   !    do y=1,Ny; y_plu=modulo(y,Ny)+1; z_plu=modulo(z,Nz)+1
   !    do x=1,Nx; x_plu=modulo(x,Nx)+1
   !       s_dp=s_dp-su2_realtrace(su2_mult_oodd(U(:,1,x,y,z),U(:,2,x_plu,y,z),U(:,1,x,y_plu,z),U(:,2,x,y,z)))
   !    end do ! x=1:Nx
   !    end do ! y=1:Ny
   !    s_slice(z) = real(s_dp) - real(s_dp_old)
   !    s_dp_old = s_dp
   !    end do ! z=1:Nz
   !    !$OMP END PARALLEL DO
   !    do z=0,Nz-1
   !       corr(z) = sum(s_slice*cshift(s_slice,z))/Nz
   !    end do
   !    s_dp = dble(1)+s_dp/dble(2*Nx*Ny*Nz)
   ! end subroutine calc_swil_corr_2D

   ! subroutine online_avgvar_dp(n_new,x,avg,var)
   !    implicit none
   !    integer,intent(in) :: n_new
   !    real(kind=sp),intent(in) :: x
   !    real(kind=dp),intent(inout) :: avg,var
   !    real(kind=dp) :: delta,obj
   !    !!! note: avg_new=((n-1)*avg_old+x)/n=avg_old+(x-avg_old)/n
   !    delta=real(x,kind=dp)-avg
   !    avg=avg+delta/dble(n_new)
   !    obj=real(x,kind=dp)-avg
   !    var=var+delta*obj !!! note: (x-avg_old)*(x-avg_new)=delta*obj
   !    !!! note: var must be divided by n_meas-1 in the end
   ! end subroutine online_avgvar_dp

#ifdef never
   pure function myavg_sp(x)
      implicit none
      real(kind=sp),dimension(:),intent(in) :: x
      real(kind=sp) :: myavg_sp,preli_sp
      integer :: length
      length=size(x,dim=1)
      preli_sp=sum(x(:))/float(length) !!! preliminary average
      myavg_sp=sum(x(:)-preli_sp)/float(length)+preli_sp
   end function myavg_sp
   pure function myavg_dp(x)
      implicit none
      real(kind=dp),dimension(:),intent(in) :: x
      real(kind=dp) :: myavg_dp,preli_dp
      integer :: length
      length=size(x,dim=1)
      preli_dp=sum(x(:))/dble(length) !!! preliminary average
      myavg_dp=sum(x(:)-preli_dp)/dble(length)+preli_dp
   end function myavg_dp

   pure function mysum_sp(x)
      implicit none
      real(kind=sp),dimension(:),intent(in) :: x
      real(kind=sp) :: mysum_sp,preli_sp
      integer :: length
      length=size(x,dim=1)
      preli_sp=sum(x(:))/float(length) !!! preliminary average
      mysum_sp=sum(x(:)-preli_sp)+preli_sp*float(length)
   end function mysum_sp
   pure function mysum_dp(x)
      implicit none
      real(kind=dp),dimension(:),intent(in) :: x
      real(kind=dp) :: mysum_dp,preli_dp
      integer :: length
      length=size(x,dim=1)
      preli_dp=sum(x(:))/dble(length) !!! preliminary average
      mysum_dp=sum(x(:)-preli_dp)+preli_dp*dble(length)
   end function mysum_dp

   pure function myavg_power_sp(x,p)
      implicit none
      real(kind=sp),dimension(:),intent(in) :: x
      integer,intent(in) :: p
      real(kind=sp) :: myavg_power_sp
      real(kind=sp),dimension(:),allocatable :: xpower
      allocate(xpower,mold=x)
      xpower(:)=x(:)**p
      myavg_power_sp=myavg_sp(xpower)
      deallocate(xpower)
   end function myavg_power_sp
   pure function myavg_power_dp(x,p)
      implicit none
      real(kind=dp),dimension(:),intent(in) :: x
      integer,intent(in) :: p
      real(kind=dp) :: myavg_power_dp
      real(kind=dp),dimension(:),allocatable :: xpower
      allocate(xpower,mold=x)
      xpower(:)=x(:)**p
      myavg_power_dp=myavg_dp(xpower)
      deallocate(xpower)
   end function myavg_power_dp

   pure function kahantwosum_sp(x)
      !!! https://accurate-algorithms.readthedocs.io/en/latest/ch04summation.html
      implicit none
      real(kind=sp),dimension(:),intent(in) :: x
      real(kind=sp) :: kahantwosum_sp,y,s,e
      integer :: lower,upper,i
      lower=lbound(x,dim=1)
      upper=ubound(x,dim=1)
      kahantwosum_sp=0.0; e=0.0
      do i=lower,upper
         y=x(i)+e !!! note: compensation step
         call twosum_sp(kahantwosum_sp,y,s,e); kahantwosum_sp=s
      end do
   end function kahantwosum_sp
   pure subroutine twosum_sp(a,b,s,e)
      !!! https://accurate-algorithms.readthedocs.io/en/latest/ch04summation.html
      implicit none
      real(kind=sp),intent( in) :: a,b
      real(kind=sp),intent(out) :: s,e
      real(kind=sp) :: z,t
      s=a+b; z=s-a
      !!! next: e=(a-(s-z))+(b-z)
      t=s-z; e=a-t
      t=b-z; e=e+t
   end subroutine twosum_sp
   pure function kahansum_sp(x)
      implicit none
      real(kind=sp),dimension(:),intent(in) :: x
      real(kind=sp) :: kahansum_sp,y,t,c,d
      integer :: lower,upper,i
      lower=lbound(x,dim=1)
      upper=ubound(x,dim=1)
      kahansum_sp=0.0; c=0.0
      do i=lower,upper
         y=x(i)-c
         t=kahansum_sp+y
         d=t-kahansum_sp
         c=d-y
         kahansum_sp=t
      end do
   end function kahansum_sp

   function jackerror_mean_sp(x,binwidth)
      implicit none !!! note: calculates jackknife-error of "mean"
      real(kind=sp),dimension(:),intent(in) :: x
      integer,intent(in) :: binwidth
      real(kind=sp) :: jackerror_mean_sp
      real(kind=sp),dimension(:),allocatable :: x_red,y
      real(kind=sp) :: x_err,y_avg,y_var
      integer :: num_bins,num_meas,num_skip,n
      num_meas=size(x,dim=1)
      num_bins=num_meas/binwidth
      num_skip=num_meas-num_bins*binwidth
      allocate(x_red((num_bins-1)*binwidth),y(num_bins))
      x_red(:)=x(num_skip+binwidth+1:num_meas) !!! note: all data but the first block and the skipped ones
      y(1)=myavg_sp(x_red(:))    !!! avg of 1st block complement
      do n=2,num_bins
         x_red((n-2)*binwidth+1:(n-1)*binwidth)=x(num_skip+(n-2)*binwidth+1:num_skip+(n-1)*binwidth) !!! note: overwrite relevant block
         y(n)=myavg_sp(x_red(:)) !!! avg of nth block complement
      end do
      y_avg=myavg_sp(y(:))
      y(:)=(y(:)-y_avg)**2
      y_var=mysum_sp(y(:))/float(num_bins-1)
      x_err=sqrt(y_var/float(num_bins))*float(num_bins-1)
      deallocate(x_red,y)
      jackerror_mean_sp=x_err
   end function jackerror_mean_sp

   function jackerror_susc_sp(x,binwidth)
      implicit none !!! note: calculates jackknife-error of "susc"
      real(kind=sp),dimension(:),intent(in) :: x
      integer,intent(in) :: binwidth
      real(kind=sp) :: jackerror_susc_sp
      real(kind=sp),dimension(:),allocatable :: x_red,x_squ,y
      real(kind=sp) :: x_err,y_avg,y_var
      integer :: num_bins,num_meas,num_skip,n
      num_meas=size(x,dim=1)
      num_bins=num_meas/binwidth
      num_skip=num_meas-num_bins*binwidth
      allocate(x_red((num_bins-1)*binwidth),x_squ((num_bins-1)*binwidth),y(num_bins))
      x_red(:)=x(num_skip+binwidth+1:num_meas) !!! note: all data but the first block and the skipped ones
      x_squ(:)=(x_red(:)-myavg_sp(x_red(:)))**2    !!! note: do not trade for x_red^2-myavg(x_red)^2
      y(1)=myavg_sp(x_squ(:))                      !!! note: avg of 1st block complement
      do n=2,num_bins
         x_red((n-2)*binwidth+1:(n-1)*binwidth)=x(num_skip+(n-2)*binwidth+1:num_skip+(n-1)*binwidth) !!! note: overwrite relevant block
         x_squ(:)=(x_red(:)-myavg_sp(x_red(:)))**2 !!! note: do not trade for x_red^2-myavg(x_red)^2
         y(n)=myavg_sp(x_squ(:))                   !!! note: avg of nth block complement
      end do
      y_avg=myavg_sp(y(:))
      y(:)=(y(:)-y_avg)**2
      y_var=mysum_sp(y(:))/float(num_bins-1)
      x_err=sqrt(y_var/float(num_bins))*float(num_bins-1)
      deallocate(x_red,x_squ,y)
      jackerror_susc_sp=x_err
   end function jackerror_susc_sp

   function jackerror_binder_sp(x,binwidth)
      implicit none !!! note: calculates jackknife-error of "binder"
      real(kind=sp),dimension(:),intent(in) :: x
      integer,intent(in) :: binwidth
      real(kind=sp) :: jackerror_binder_sp
      real(kind=sp),dimension(:),allocatable :: x_red,x_squ,x_fou,y
      real(kind=sp) :: x_err,y_avg,y_var
      integer :: num_bins,num_meas,num_skip,n
      num_meas=size(x,dim=1)
      num_bins=num_meas/binwidth
      num_skip=num_meas-num_bins*binwidth
      allocate(x_red((num_bins-1)*binwidth),x_squ((num_bins-1)*binwidth),x_fou((num_bins-1)*binwidth),y(num_bins))
      x_red(:)=x(num_skip+binwidth+1:num_meas) !!! note: all data but the first block and the skipped ones
      x_squ(:)=x_red(:)**2
      x_fou(:)=x_red(:)**4
      y(1)=1.0-myavg_sp(x_fou(:))/myavg_sp(x_squ)**2/3.0    !!! note: avg of 1st block complement
      do n=2,num_bins
         x_red((n-2)*binwidth+1:(n-1)*binwidth)=x(num_skip+(n-2)*binwidth+1:num_skip+(n-1)*binwidth) !!! note: overwrite relevant block
         x_squ(:)=x_red(:)**2
         x_fou(:)=x_red(:)**4
         y(n)=1.0-myavg_sp(x_fou(:))/myavg_sp(x_squ)**2/3.0 !!! note: avg of nth block complement
      end do
      y_avg=myavg_sp(y(:))
      y(:)=(y(:)-y_avg)**2
      y_var=mysum_sp(y(:))/float(num_bins-1)
      x_err=sqrt(y_var/float(num_bins))*float(num_bins-1)
      deallocate(x_red,x_squ,x_fou,y)
      jackerror_binder_sp=x_err
   end function jackerror_binder_sp
#endif

!   pure elemental function psqu_naive(x,y)
!      implicit none
!      integer,intent(in) :: x,y
!      real(kind=sp) :: psqu_naive
!      psqu_naive=(twopi_sp*float(x)/float(Nx))**2+(twopi_sp*float(y)/float(Ny))**2
!   end function psqu_naive
!   pure elemental function psqu_boson(x,y)
!      implicit none
!      integer,intent(in) :: x,y
!      real(kind=sp) :: psqu_boson
!      psqu_boson=4.0*sin(0.5*twopi_sp*float(x)/float(Nx))**2+4.0*sin(0.5*twopi_sp*float(y)/float(Ny))**2
!   end function psqu_boson
!   pure function find_firstnegative(obj)
!      implicit none
!      real(kind=dp),intent(in) :: obj(0:Ny/2)
!      integer :: find_firstnegative,y
!      find_firstnegative=huge(find_firstnegative)
!      do y=Ny/2,0,-1
!         if (obj(y)<0.0_dp) find_firstnegative=y
!      end do
!   end function find_firstnegative
!   pure function effmass_naive(corr_tmin,corr_tmaj)
!      implicit none
!      real(kind=dp),intent(in) :: corr_tmin,corr_tmaj
!      real(kind=sp) :: effmass_naive
!      if ((corr_tmin>0.0).and.(corr_tmaj>0.0).and.(corr_tmin/corr_tmaj>1.0)) then
!         effmass_naive=0.5*sngl(max(log(corr_tmin/corr_tmaj),0.0_dp))
!      else
!         effmass_naive=0.0
!      end if
!   end function effmass_naive
!   pure function effmass_robust(corr_tmin,corr_t,corr_tmaj)
!      implicit none
!      real(kind=dp),intent(in) :: corr_tmin,corr_t,corr_tmaj
!      real(kind=sp) :: effmass_robust
!      if (corr_t>0.0) then
!         effmass_robust=sngl(sqrt(max((corr_tmin-2.0_dp*corr_t+corr_tmaj)/(corr_t),0.0_dp)))
!      else
!         effmass_robust=0.0
!      end if
!   end function effmass_robust
!   pure function effmass_balint(corr_tmin,corr_tmaj,corr_midpoint)
!      implicit none
!      real(kind=dp),intent(in) :: corr_tmin,corr_tmaj,corr_midpoint
!      real(kind=sp) :: effmass_balint
!      if ((corr_tmin<0.0).or.(corr_tmaj<0.0).or.(corr_midpoint<0.0)) then
!         effmass_balint=0.0
!      else if ((corr_tmin/corr_midpoint>1.0).and.(corr_tmaj/corr_midpoint>1.0)) then
!         effmass_balint=0.5*sngl(max(acosh(corr_tmin/corr_midpoint)-acosh(corr_tmaj/corr_midpoint),0.0_dp))
!      else
!         effmass_balint=0.0
!      end if
!   end function effmass_balint

end program su2_3D
