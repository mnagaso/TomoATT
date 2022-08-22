include "eikon_solver_mpi.f90"


! example: 3-D 3rd order anisotropic eikonal equation in cartesian coordinate (Point source)
! parameter setting:
! domain:                   R * Theta * Phi = [6275, 6375] * [49^\circ, 51^\circ] * [129^\circ, 131^\circ]
! analytic solution:        T = |x-x_s|/c0
! isotropic eik equ:      Tx^2 + Ty^2 + Tz^2 = s^2,
! boundary condition:       T(x0,y0,z0) = 0  (point source)
! test function:    'FSM_WENO3_PS_sphe_3d' in "eikon_solver.f90"
!

program eikonal_2d
    use mpi

! ######################### 参数声明 parameters statement #########################

    implicit none

    ! mesh grid
    integer,parameter :: nr=55,nt=55,np=55
    double precision,parameter :: pi=3.14159265358979323846264338327950288
    double precision,parameter :: rr1=6070,rr2=6400
    double precision,parameter :: tt1=(30.0-1.5)/180*pi,tt2=(50.0+1.5)/180*pi
    double precision,parameter :: pp1=(15.0-1.5)/180*pi,pp2=(40.0+1.5)/180*pi
    double precision :: rr(nr),tt(nt),pp(np),dr,dt,dp


    ! source (STATIONS) & receiver (EARTHQUAKES)
    integer :: nsou
    double precision,allocatable :: rsou(:),tsou(:),psou(:)  ! source (station)
    integer :: nrec
    double precision,allocatable :: rrec(:),trec(:),prec(:)    ! receiver (earthquake)

    ! model parameter
    double precision,parameter :: gamma = 0.0
    double precision :: xi(nr,nt,np),eta(nr,nt,np),zeta(nr,nt,np)   ! syn model (ani)
    double precision :: xiT(nr,nt,np),etaT(nr,nt,np),zetaT(nr,nt,np)   ! true model (ani)
    double precision :: sigma1,sigma2,psi
    double precision,parameter :: slow_p=0.04,ani_p=0.03

    ! Eikonal solver
    double precision :: aT(nr,nt,np),bT(nr,nt,np),cT(nr,nt,np),fT(nr,nt,np),funT(nr,nt,np)  ! true model 真实模型
    double precision :: a(nr,nt,np),b(nr,nt,np),c(nr,nt,np),f(nr,nt,np),fun(nr,nt,np)   ! syn model 猜测模型
    double precision :: uT(nr,nt,np),TableT(nr,nt,np)       !true timetable 真实走时场
    double precision :: u(nr,nt,np),Table(nr,nt,np)         !syn timetable 猜测走时场
    !double precision :: T0para(8),tau(nr,nt,np)    ! multiplicative factor

    ! adjoint source
    double precision,allocatable :: sourceADJ(:),TtimeT(:,:),Ttime(:,:) ! adjoint source 伴随源,真实走时，合成走时
    double precision :: TableADJ(nr,nt,np)  ! adjoint field 伴随场
    double precision :: Ks(nr,nt,np),Kxi(nr,nt,np),Keta(nr,nt,np)  ! event kernel
    integer,parameter :: nk = 3   ! number of kernels
    double precision :: all_Ks(nr,nt,np),all_Kxi(nr,nt,np),all_Keta(nr,nt,np)  ! sensitivity kernel
    double precision :: all_Kernel(nr,nt,np,nk) ! kernels for all parameters together
    !double precision,parameter :: radius = 40 ! kernel mask radius

    ! inversion grid
    integer,parameter :: ninvr=11,ninvt=10,ninvp=10   ! inversion grid size (coarse)
    double precision,parameter :: invr1=6070.0,invr2=6400.0
    double precision,parameter :: invt1=30.0/180*pi,invt2=50.0/180*pi
    double precision,parameter :: invp1=15.0/180*pi,invp2=40.0/180*pi
    integer,parameter :: ngrid=5  ! number of nultigrid (multigrid technique)
    double precision :: invr(ngrid,ninvr),invt(ngrid,ninvt),invp(ngrid,ninvp),dinvr,dinvt,dinvp
    double precision :: update_value(nr,nt,np,nk)   ! parameter update value 更新参数的变化量

    ! model update (optimization)
    double precision :: stepsize    ! stepsize
    integer,parameter :: MaxIter=1  ! max iteration step
    double precision :: old_obj, obj, sou_obj ! 优化目标函数 objective function
    double precision,allocatable :: weight_sou(:),weight_rec(:) !weight function
    double precision :: dis_ij, dis_zero     !weight function
    logical :: weight_sou_bool=.false. , weight_rec_bool= .false.

    ! temp var
    double precision :: depth,lat,lon

    ! output file
    character(Len=80) :: filename,form,form2

    ! loop index
    integer :: iir,iit,iip,idi,idj,idk,i,j,k,iter,igrid

    ! computing time
    double precision :: time_begin,time_end,alpha

    ! mpi parameter
    integer :: ierr,myid,nproc

    ! #############################  main work #############################

    call mpi_init(ierr)

    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    call CPU_TIME(time_begin)   ! program start



! ################### read sources and receivers 读取地震台站数据 ###################
    open(10,file='ega5/input/STATIONS_5.txt')
    read(10,*) nsou
    allocate(rsou(nsou),tsou(nsou),psou(nsou),weight_sou(nsou))
    do i=1,nsou
        read(10,*) depth,lat,lon
        rsou(i) = 6371.0-depth
        tsou(i) = lat/180.0*pi
        psou(i) = lon/180.0*pi
    end do
    close(10)

    open(10,file='ega5/input/EARTHQUAKES_5.txt')
    read(10,*) nrec
    allocate(rrec(nrec),trec(nrec),prec(nrec),weight_rec(nrec))
    do i=1,nrec
        read(10,*) depth,lat,lon
        rrec(i) = 6371.0-depth
        trec(i) = lat/180.0*pi
        prec(i) = lon/180.0*pi
    end do
    close(10)

    ! distance based weight -rec 台站权重
    if (weight_rec_bool) then
        dis_zero = 0
        do i=1,nrec-1
            do j=i,nrec
                call Epicenter_Distance(trec(i),prec(i),trec(j),prec(j),dis_ij)
                dis_zero = dis_zero + dis_ij
            end do
        end do
        dis_zero = dis_zero/nrec/(nrec-1)*2
        dis_zero = dis_zero/2

        do i=1,nrec
            weight_rec(i)=0.0
            do j=1,nrec
                call Epicenter_Distance(trec(i),prec(i),trec(j),prec(j),dis_ij)
                weight_rec(i) = weight_rec(i) + exp(-(dis_ij/dis_zero)**2)
            end do
            weight_rec(i) = 1.0/weight_rec(i)
        end do
        print *, dis_zero
    end if

    if (weight_rec_bool) then
    ! distance based weight - soruce 地震权重
        dis_zero = 0
        do i=1,nsou-1
            do j=i,nsou
                call Epicenter_Distance(tsou(i),psou(i),tsou(j),psou(j),dis_ij)
                dis_zero = dis_zero + dis_ij
            end do
        end do
        dis_zero = dis_zero/nsou/(nsou-1)*2
        dis_zero = dis_zero/4

        do i=1,nsou
            weight_sou(i)=0.0
            do j=1,nsou
                call Epicenter_Distance(tsou(i),psou(i),tsou(j),psou(j),dis_ij)
                weight_sou(i) = weight_sou(i) + exp(-(dis_ij/dis_zero)**2)
            end do
            weight_sou(i) = 1.0/weight_sou(i)
        end do
        print *, dis_zero
    end if

    ! ---- 记录到 matlab 数据中 matlab record ----
    if (myid .eq. 0) then
        open(100,file='ega5/output/STATIONS')
        open(101,file='ega5/output/EARTHQUAKES')
            do i=1,nsou
                write(100,'(f13.7,f13.7,f13.7,f13.7)') rsou(i),tsou(i),psou(i),weight_sou(i)
            end do
            do i=1,nrec
                write(101,'(f13.7,f13.7,f13.7,f13.7)') rrec(i),trec(i),prec(i),weight_rec(i)
            end do
        close(100);close(101)
    end if

! #########################  生成网格 build mesh grid  #####################################
    ! -------- forward modeling grid ----------
    dr=(rr2-rr1)/(nr-1); dt=(tt2-tt1)/(nt-1); dp=(pp2-pp1)/(np-1)
    do iir=1,nr
        rr(iir)=rr1+(iir-1)*dr
    end do
    do iit=1,nt
        tt(iit)=tt1+(iit-1)*dt
    end do
    do iip=1,np
        pp(iip)=pp1+(iip-1)*dp
    end do

    ! -------- inversion multigrid ----------
    dinvr=(invr2-invr1)/(ninvr-1); dinvt=(invt2-invt1)/(ninvt-1); dinvp=(invp2-invp1)/(ninvp-1)
    do igrid=1,ngrid
        do iir=1,ninvr
            invr(igrid,iir)=invr1+(iir-1)*dinvr-(igrid-1)*dinvr/ngrid
        end do
        do iit=1,ninvt
            invt(igrid,iit)=invt1+(iit-1)*dinvt-(igrid-1)*dinvt/ngrid
        end do
        do iip=1,ninvp
            invp(igrid,iip)=invp1+(iip-1)*dinvp-(igrid-1)*dinvp/ngrid
        end do
    end do

! #########################  构建背景模型 build the synthetic and true model #####################################

    do iir=1,nr
        do iit=1,nt
            do iip=1,np

                ! synthetic (initial) model
                eta(iir,iit,iip)=0.0;
                xi(iir,iit,iip)=0.0;
                zeta(iir,iit,iip)=gamma*sqrt(eta(iir,iit,iip)**2+xi(iir,iit,iip)**2)
                ! if (rr(iir)>6351) then      ! 6371 - 6351   20 km
                !     fun(iir,iit,iip) = 1.0/(5.8+(6371-rr(iir))/20.0*0.7)
                ! elseif (rr(iir)>6336) then  ! 6351 - 6336   15 km
                !     fun(iir,iit,iip) = 1.0/(6.5+(6351-rr(iir))/15.0*0.6)
                ! elseif (rr(iir)>6251) then  ! 6336 - 6251   85 km
                !     fun(iir,iit,iip) = 1.0/(8.04+(6336-rr(iir))/85.0*0.01)
                ! elseif (rr(iir)>6161) then  ! 6251 - 6161   90 km
                !     fun(iir,iit,iip) = 1.0/(8.05+(6251-rr(iir))/90.0*0.25)
                ! elseif (rr(iir)>5961) then  ! 6161 - 5961   200 km
                !     fun(iir,iit,iip) = 1.0/(8.30+(6161-rr(iir))/200.0*0.73)
                ! else
                !     fun(iir,iit,iip) = 1.0/9.03
                ! end if

                if (rr(iir)>6351) then      ! 6371 - 6351   20 km
                    fun(iir,iit,iip) = 1.0/(5.8+(6371-rr(iir))/20.0*0.7)
                elseif (rr(iir)>6336) then  ! 6351 - 6336   15 km
                    fun(iir,iit,iip) = 1.0/(6.5+(6351-rr(iir))/15.0*0.6)
                elseif (rr(iir)>5961) then  ! 6351 - 6336   15 km
                    fun(iir,iit,iip) = 1.0/(8.0+(6336-rr(iir))/375.0*1)
                else
                    fun(iir,iit,iip) = 1.0/9.0
                end if


                !read(100,*) xi(iir,iit,iip)
                !read(101,*) eta(iir,iit,iip)
                !read(102,*) fun(iir,iit,iip)

                a(iir,iit,iip)=1.0+2*zeta(iir,iit,iip);
                b(iir,iit,iip)=1.0-2*xi(iir,iit,iip);
                c(iir,iit,iip)=1.0+2*xi(iir,iit,iip);
                f(iir,iit,iip)=-2*eta(iir,iit,iip);

                Table(iir,iit,iip) = 0



                ! true (target) model
                if (tt(iit)>=30.0/180*pi .and. tt(iit)<=50.0/180*pi .and. &
                  & pp(iip)>=15.0/180*pi .and. pp(iip)<=40.0/180*pi .and. &
                  & rr(iir)>=6211.0      .and. rr(iir)<=6371.0 ) then
                    sigma1 = sin(4*pi*(tt(iit)-30.0/180*pi)/(20.0/180*pi))* &
                           & sin(4*pi*(pp(iip)-15.0/180*pi)/(25.0/180*pi))* &
                           & sin(2*pi*(rr(iir)-6211)/160.0)
                else
                    sigma1 = 0.0
                end if


                if (sigma1<0) then
                    psi = 60.0/180*pi
                elseif (sigma1 > 0) then
                    psi = 150.0/180*pi
                else
                    psi = 0
                end if

                etaT(iir,iit,iip)=ani_p*abs(sigma1)*sin(2*psi);
                xiT(iir,iit,iip)=ani_p*abs(sigma1)*cos(2*psi);
                zetaT(iir,iit,iip)=gamma*sqrt(etaT(iir,iit,iip)**2+xiT(iir,iit,iip)**2)

                aT(iir,iit,iip)=1.0+2*zetaT(iir,iit,iip);
                bT(iir,iit,iip)=1.0-2*xiT(iir,iit,iip);
                cT(iir,iit,iip)=1.0+2*xiT(iir,iit,iip);
                fT(iir,iit,iip)=-2*etaT(iir,iit,iip);
                funT(iir,iit,iip) = fun(iir,iit,iip)/(1+sigma1*slow_p)

                !fun(iir,iit,iip) = funT(iir,iit,iip)
                !eta(iir,iit,iip)=etaT(iir,iit,iip);
                !xi(iir,iit,iip)=xiT(iir,iit,iip);
                !zeta(iir,iit,iip)=zetaT(iir,iit,iip)
                !a(iir,iit,iip)=aT(iir,iit,iip)
                !b(iir,iit,iip)=bT(iir,iit,iip)
                !c(iir,iit,iip)=cT(iir,iit,iip)
                !f(iir,iit,iip)=fT(iir,iit,iip)

                TableT(iir,iit,iip) = 0
                u(iir,iit,iip) = 0

            end do
        end do
    end do



    ! ---- 记录到 matlab 数据中 matlab record ----
    if (myid==0) then
        open(100,file='ega5/output/model_true')
        open(101,file='ega5/output/model_init')
            do iir=1,nr
                do iit=1,nt
                    do iip=1,np
                        write(100,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)')  &
                            & rr(iir),tt(iit),pp(iip), &
                            & funT(iir,iit,iip),xiT(iir,iit,iip),etaT(iir,iit,iip)
                        write(101,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)')  &
                            & rr(iir),tt(iit),pp(iip), &
                            & fun(iir,iit,iip),xi(iir,iit,iip),eta(iir,iit,iip)
                    end do
                end do
            end do
        close(100);close(101)

        ! --------- 多重网格 multiple grid ------
        do igrid=1,Ngrid
            select case(igrid)
            case (1:9)
                write(form,'(I1)') igrid
            case (10:99)
                write(form,'(I2)') igrid
            end select
            open(100,file='ega5/output/multigrid'//trim(form))

            do iir=1,ninvr
                do iit=1,ninvt
                    do iip=1,ninvp
                        write(100,'(f13.7,f13.7,f13.7)') &
                        & invr(igrid,iir),invt(igrid,iit),invp(igrid,iip)
                    end do
                end do
            end do
            close(100)
        end do

    end if


! ######################### 计算真实到时 calculate true traveltimes ######################
    call mpi_barrier(mpi_comm_world,ierr)
    allocate(TtimeT(nsou,nrec)) ! 初始化真实到时表 initiate true traveltimes
    allocate(Ttime(nsou,nrec))  ! 初始化合成到时表 initiate syn traveltimes
    allocate(sourceADJ(nrec))   ! 初始化伴随源  initiate adjoint source

    if (myid==0) then
        print *, ' '
        print *, '----------------- calculating true timetable ... ----------------'
        print *, ' '
    end if

    do i=1,nsou  ! loop sources (STATIONS)
        ! --------  求解程函方程 solve eikonal equations ------------
        call mpi_barrier(mpi_comm_world,ierr)
        !rsou(i)=6371.0; tsou(i)=30.5/180.0*pi; psou(i)=15.5/180.0*pi
        call FSM_WENO3_PS_sphe_3d_mul_mpi(rr,tt,pp,nr,nt,np,aT,bT,cT,fT,TableT,funT,rsou(i),tsou(i),psou(i),u)
        call mpi_barrier(mpi_comm_world,ierr)
        ! --------  得到真实走时 calculate true arrival time at receiver --------
        call Linear_Interp_3D(rr,tt,pp,TableT,nr,nt,np,rrec,trec,prec,TtimeT(i,:),nrec)
        !if ((myid.eq.0) .and. 1>0) then
        !    print *, rsou(i),tsou(i)/pi*180,psou(i)/pi*180
        !    open(100,file='ega5/output/arrivaltime')
        !    do j=1,nrec
        !        write(100,'(i5.3,f8.2,f6.2,f6.2,f8.3)') j, rrec(j),trec(j)/pi*180,prec(j)/pi*180,TtimeT(i,j)
        !    end do
        !    close(100)
        !end if
    end do

    print *, '----------------- calculating true timetable over ---------------- '

! ########################  反演开始  inversion start ##########################

    ! ----- 初始化反演参数 initiate inversion parameters ----
    stepsize = 0.01 ! 迭代步长 stepsize
    old_obj = 0; obj = 0 ! 目标函数 objective function



    if (myid .eq. 0) then
        print *, ' '
        print *, '----------------- inversion start ... ----------------'
        print *, ' '
    end if



    do iter = 1,MaxIter
        call mpi_barrier(mpi_comm_world,ierr)
        if (myid .eq. 0) then
            print *, '----------------- iteration ',iter,' starting ... ----------------'
        end if

        ! ----- 初始化参数 initiate parameters ------
        all_Ks = 0; all_Kxi = 0; all_Keta = 0;
        obj = 0;

        if (myid .eq. 0) then
            select case (iter)
            case (1:9)
                write(form,'(I1)') iter
            case (10:99)
                write(form,'(I2)') iter
            end select
            filename='ega5/output/misfit_step'//trim(form)
            open(998, file=trim(filename))
        end if


        do i=1,nsou  ! loop sources (STATIONS)

! ########################  计算合成走时场 calculate synthetic timetable ########################
            if (myid .eq. 0) then
                print *, '----------------- calculating synthetic timetable ... ----------------'
            end if

            Table = 20
            call mpi_barrier(mpi_comm_world,ierr)
            call FSM_WENO3_PS_sphe_3d_mul_mpi(rr,tt,pp,nr,nt,np,a,b,c,f,Table,fun,rsou(i),tsou(i),psou(i),u)
            call mpi_barrier(mpi_comm_world,ierr)
            ! --------  得到合成走时 calculate synthetic arrival time at receiver --------
            call Linear_Interp_3D(rr,tt,pp,Table,nr,nt,np,rrec,trec,prec,Ttime(i,:),nrec)

            ! --------  构造伴随源 build adjoint source -------
            call Adjoint_Source_Dt(Ttime(i,:),TtimeT(i,:),nrec,sourceADJ,sou_obj)  ! absolute traveltime difference
            if (myid .eq. 0) then
                do j=1,nrec
                    write(998,'(f10.5)') Ttime(i,j)-TtimeT(i,j)
                end do
            end if

            ! weighting adjoint source
            if (weight_rec_bool) then
                do j=1,nrec
                    sourceADJ(j) = sourceADJ(j) * weight_rec(j)
                end do
            end if

            if (1>0 .and. (myid .eq. 0)) then
                select case (iter)
                case (1:9)
                    write(form,'(I1)') iter
                case (10:99)
                    write(form,'(I2)') iter
                end select
                select case (i)
                case (1:9)
                    write(form2,'(I1)') i
                case (10:99)
                    write(form2,'(I2)') i
                end select
                filename='ega5/output/syn_step'//trim(form)//'_event'//trim(form2)
                open(100,file=trim(filename))
                do iir=1,nr
                    do iit=1,nt
                        do iip=1,np
                            write(100,'(f13.7)') Table(iir,iit,iip)
                        end do
                    end do
                end do
                close(100)
            end if
!

            !print *, Ttime(i,:),TtimeT(i,:)
            !call Adjoint_Source_Ddt(Ttime(i,:),TtimeT(i,:),nrec,sourceADJ,sou_obj)  ! double difference traveltime
            !do j=1,nrec-1
            !    do k=j,nrec
            !        write(998,'(f10.5)') (Ttime(i,j)-Ttime(i,k))-(TtimeT(i,j)-TtimeT(i,k))
            !   end do
            !end do

! ######################## 计算伴随场 adjoint field ########################
            if (myid .eq. 0) then
            !    print *, '----------------- calculating adjoint field ... ----------------'
            end if
            call FSM_O1_Adj_sphe_3d(rr,tt,pp,nr,nt,np,Table,TableADJ,zeta,xi,eta,nrec,rrec,trec,prec,sourceADJ)
            !call FSM_O1_Adj_PS_sphe_3d(rr,tt,pp,nr,nt,np,T0para,tau,TableADJ,zeta,xi,eta,nrec,rrec,trec,prec,sourceADJ)
            ! ---- 记录到 matlab 数据中 matlab record (adjoint field) ----
            if (1>0 .and. (myid .eq. 0)) then
                select case (iter)
                case (1:9)
                    write(form,'(I1)') iter
                case (10:99)
                    write(form,'(I2)') iter
                end select
                select case (i)
                case (1:9)
                    write(form2,'(I1)') i
                case (10:99)
                    write(form2,'(I2)') i
                end select
                filename='ega5/output/adj_step'//trim(form)//'_event'//trim(form2)
                open(100,file=trim(filename))
                do iir=1,nr
                    do iit=1,nt
                        do iip=1,np
                            write(100,'(f13.7)') TableADJ(iir,iit,iip)
                        end do
                    end do
                end do
                close(100)
            end if
! ######################## 计算敏感核 sensitivity kernel ########################
            !print *, '----------------- calculating kernel ... ----------------'
            call Sensitivity_Kernel(rr,tt,pp,nr,nt,np,Table,TableADJ,gamma,xi,eta,fun,Ks,Kxi,Keta)

            ! ------- 抹去伴随源处的大函数值 mask the source -------
            !print *, '----------------- kernel mask ----------------'
            !call Kernel_Mask(rr,tt,pp,nr,nt,np,Ks,rsou(i),tsou(i),psou(i),radius)
            !call Kernel_Mask(rr,tt,pp,nr,nt,np,Kxi,rsou(i),tsou(i),psou(i),radius)
            !call Kernel_Mask(rr,tt,pp,nr,nt,np,Keta,rsou(i),tsou(i),psou(i),radius)
            call Kernel_Mask_new(rr,tt,pp,nr,nt,np,Ks,rsou(i),tsou(i),psou(i))
            call Kernel_Mask_new(rr,tt,pp,nr,nt,np,Kxi,rsou(i),tsou(i),psou(i))
            call Kernel_Mask_new(rr,tt,pp,nr,nt,np,Keta,rsou(i),tsou(i),psou(i))

            ! ---- 记录到 matlab 数据中 matlab record (adjoint field) ----
            if (1<0 .and. myid .eq. 0) then
                select case (iter)
                case (1:9)
                    write(form,'(I1)') iter
                case (10:99)
                    write(form,'(I2)') iter
                end select
                select case (i)
                case (1:9)
                    write(form2,'(I1)') i
                case (10:99)
                    write(form2,'(I2)') i
                end select
                filename='ega5/output/kernel_step'//trim(form)//'_event'//trim(form2)
                open(100,file=trim(filename))

                do iir=1,nr
                    do iit=1,nt
                        do iip=1,np
                            write(100,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)') &
                            & rr(iir),tt(iit),pp(iip), &
                            & Ks(iir,iit,iip),Kxi(iir,iit,iip),Keta(iir,iit,iip)
                        end do
                    end do
                end do
                close(100)
            end if

            ! --------- 敏感核,obj叠加 sum kernels and objective function -------
            if (weight_sou_bool) then   ! no weight
                all_Ks = all_Ks + Ks * weight_sou(i);
                all_Kxi = all_Kxi + Kxi * weight_sou(i);
                all_Keta = all_Keta + Keta * weight_sou(i);
                obj = obj + sou_obj !* weight_sou(i)
            else    ! source weighted
                all_Ks = all_Ks + Ks;
                all_Kxi = all_Kxi + Kxi;
                all_Keta = all_Keta + Keta;
                obj = obj + sou_obj
            end if

        end do

        close(998)

        ! ---- 记录到 matlab 数据中 matlab record (kernel) ----
        if (1>0 .and. myid .eq. 0) then
            select case (iter)
            case (1:9)
                write(form,'(I1)') iter
            case (10:99)
                write(form,'(I2)') iter
            end select
            filename='ega5/output/kernel_step'//trim(form)//'_sum'
            open(100,file=trim(filename))
            do iir=1,nr
                do iit=1,nt
                    do iip=1,np
                        write(100,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)') &
                            & rr(iir),tt(iit),pp(iip), &
                            & all_Ks(iir,iit,iip),all_Kxi(iir,iit,iip),all_Keta(iir,iit,iip)
                    end do
                end do
            end do
            close(100)
        end if
! ##################### 模型更新 model update ####################
        print *, '----------------- model updating ... ----------------'

        all_Kernel(:,:,:,1) = all_Ks;
        all_Kernel(:,:,:,2) = all_Kxi;
        all_Kernel(:,:,:,3) = all_Keta;

        ! ------------ 更新速度(慢度) update velocity (slowness) ------------
        call Parameter_Update_Multigrid(rr,tt,pp,nr,nt,np,all_Kernel,nk, &
                                        & invr,invt,invp,ninvr,ninvt,ninvp,ngrid,stepsize,update_value)

        if (1<0 .and. myid .eq. 0) then
            select case (iter)
            case (1:9)
                write(form,'(I1)') iter
            case (10:99)
                write(form,'(I2)') iter
            end select
            filename='ega5/output/model_update_step'//trim(form)
            open(100,file=trim(filename))
            do iir=1,nr
                do iit=1,nt
                    do iip=1,np
                        write(100,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)') &
                            & rr(iir),tt(iit),pp(iip), &
                            & -fun(iir,iit,iip)*update_value(iir,iit,iip,1), &
                            & -update_value(iir,iit,iip,2),-update_value(iir,iit,iip,3)
                    end do
                end do
            end do
            close(100)
        end if

        do iir=1,nr
            do iit=1,nt
                do iip=1,np
                    fun(iir,iit,iip) = fun(iir,iit,iip)*(1-update_value(iir,iit,iip,1))
                    !fun(iir,iit,iip) = fun(iir,iit,iip)
                end do
            end do
        end do

        ! ------------ 更新xi update xi ------------
        do iir=1,nr
            do iit=1,nt
                do iip=1,np
                    xi(iir,iit,iip) = xi(iir,iit,iip)-update_value(iir,iit,iip,2)
                end do
            end do
        end do

        ! ------------ 更新eta update eta ------------
        do iir=1,nr
            do iit=1,nt
                do iip=1,np
                    eta(iir,iit,iip) = eta(iir,iit,iip)-update_value(iir,iit,iip,3)
                end do
            end do
        end do

        ! ----------- 更新 update eikonal solver paramater -----------------
        do iir=1,nr
            do iit=1,nt
                do iip=1,np
                    b(iir,iit,iip)=1.0-2*xi(iir,iit,iip);
                    c(iir,iit,iip)=1.0+2*xi(iir,iit,iip);
                    f(iir,iit,iip)=-2*eta(iir,iit,iip);
                end do
            end do
        end do




! ##################### 调整下降步长 modify stepsize #############
        if (myid .eq. 0) then
            print *, ' '
        end if

        if (1>0 .and. (myid .eq. 0)) then
            open(999,file='ega5/output/obj',access='append')
            write(999,'(f9.2,f9.6)') obj, stepsize
            close(999)
        end if

        if (iter == 1 ) then
            if (myid .eq. 0) then
                write(*,'(a,f9.2)') 'iter 1, obj is', obj
                write(*,'(a,f9.6)') 'iter 1, stepsize is', stepsize
            end if
        elseif (iter >= 2 .and. obj < old_obj) then
            !stepsize = min(0.01,stepsize*1.2)
            if (myid .eq. 0) then
                write(*,'(a,f9.2,a,f9.2)') 'objective function decreases, from', old_obj, ' to', obj
                write(*,'(a,f9.6)') 'new stepsize is ', stepsize
            end if
        elseif (iter >= 2 .and. obj >= old_obj) then
            !stepsize = max(0.002,stepsize*0.5)
            if (myid .eq. 0) then
                write(*,'(a,f9.2,a,f9.2)') 'objective function increases, from', old_obj, ' to', obj
                write(*,'(a,f9.6)') 'new stepsize is ', stepsize
            end if
        end if

        print *, '----------------- iteration ',iter,' over ----------------'
        print *, ' '
        print *, ' '


        old_obj = obj

        ! ---- 记录到 matlab 数据中 matlab record ----
        if (1>0 .and. (myid .eq. 0)) then
            select case (iter)
            case (1:9)
                write(form,'(I1)') iter
            case (10:99)
                write(form,'(I2)') iter
            end select
            filename='ega5/output/model_step'//trim(form)
            open(100,file=trim(filename))
            do iir=1,nr
                do iit=1,nt
                    do iip=1,np
                        write(100,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)') &
                            & rr(iir),tt(iit),pp(iip), &
                            & fun(iir,iit,iip),xi(iir,iit,iip),eta(iir,iit,iip)
                    end do
                end do
            end do
            close(100)
        end if


        if (myid .eq. 0) then
            call CPU_TIME(time_end)
            write(*,'(f10.2)') time_end - time_begin
        end if

    end do


    call mpi_finalize(ierr)

end program eikonal_2d

