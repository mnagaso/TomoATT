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
    integer :: nr,nt,np
    double precision,parameter :: pi=3.14159265358979323846264338327950288
    double precision,parameter :: Earth_Radius = 6371.0
    double precision :: rr1,rr2,tt1,tt2,pp1,pp2
    double precision,allocatable :: rr(:),tt(:),pp(:)
    double precision :: dr,dt,dp

    ! STATIONS & EARTHQUAKES
    integer :: nev
    double precision,allocatable :: rev(:),tev(:),pev(:)  ! earthquake (event)
    integer :: nst
    double precision,allocatable :: rst(:),tst(:),pst(:)    ! station
    integer :: nev_tele
    double precision,allocatable :: rev_tele(:),tev_tele(:),pev_tele(:)  ! tele seismic event 

    ! model parameter
    double precision,parameter :: gamma = 0.0
    double precision :: slow_p,ani_p
    double precision,allocatable :: xi(:,:,:),eta(:,:,:),zeta(:,:,:)   ! syn model (ani)
    double precision,allocatable :: a(:,:,:),b(:,:,:),c(:,:,:),f(:,:,:),fun(:,:,:)   ! syn model 猜测模型 
    double precision,allocatable :: xiT(:,:,:),etaT(:,:,:),zetaT(:,:,:)   ! true model (ani)
    double precision,allocatable :: aT(:,:,:),bT(:,:,:),cT(:,:,:),fT(:,:,:),funT(:,:,:)  ! true model 真实模型
    double precision :: sigma1,sigma2,psi
    
    ! teleseismic modeling
    double precision,allocatable :: bd_N(:,:,:),bd_S(:,:,:),bd_W(:,:,:),bd_E(:,:,:),bd_B(:,:,:),bd_T(:,:,:) ! boundary conditions
    double precision,allocatable :: bd_N_l(:,:,:),bd_S_l(:,:,:),bd_W_l(:,:,:),bd_E_l(:,:,:),bd_B_l(:,:,:),bd_T_l(:,:,:) ! local boundary conditions
    
    logical,allocatable :: isbd(:,:)
    double precision :: max_deg,deg 

    ! tele and regional control
    double precision :: region_ratio, tele_ratio    ! the ratio of regional kernel over tele kernel
    logical :: is_region, is_tele


    ! Eikonal solver
    double precision,allocatable :: uT(:,:,:),TableT(:,:,:)       !true timetable 真实走时场
    double precision,allocatable :: u(:,:,:),Table(:,:,:)         !syn timetable 猜测走时场
    double precision,allocatable :: TableADJ(:,:,:)  ! adjoint field 伴随场

    ! adjoint source 
    double precision,allocatable :: TtimeT(:,:),Ttime_teleT(:,:),Ttime(:,:),Ttime_tele(:,:) ! 真实走时，合成走时
    double precision,allocatable :: TtimeT_l(:,:),Ttime_teleT_l(:,:),Ttime_l(:,:),Ttime_tele_l(:,:) ! 真实走时，合成走时
    double precision,allocatable :: sourceADJ(:)    ! adjoint source 伴随源
    
    ! sensitivity kernel
    integer :: nk   ! number of kernels
    double precision,allocatable :: all_Ks_region(:,:,:),all_Kxi_region(:,:,:),all_Keta_region(:,:,:)  ! sensitivity kernel
    double precision,allocatable :: all_Ks_region_l(:,:,:),all_Kxi_region_l(:,:,:),all_Keta_region_l(:,:,:)  ! sensitivity kernel
    double precision,allocatable :: all_Ks_tele(:,:,:),all_Kxi_tele(:,:,:),all_Keta_tele(:,:,:)  ! sensitivity kernel
    double precision,allocatable :: all_Ks_tele_l(:,:,:),all_Kxi_tele_l(:,:,:),all_Keta_tele_l(:,:,:)  ! sensitivity kernel
    double precision,allocatable :: Ks(:,:,:),Kxi(:,:,:),Keta(:,:,:)  ! event kernel
    double precision,allocatable :: all_Kernel(:,:,:,:) ! kernels for all parameters together

    ! inversion parameterization
    integer :: ninvr,ninvt,ninvp
    double precision :: invr1,invr2,invt1,invt2,invp1,invp2
    integer :: ngrid  ! number of nultigrid (multigrid technique) 
    double precision,allocatable :: invr(:,:),invt(:,:),invp(:,:)
    double precision,allocatable :: inv_dep(:)
    double precision :: dinvr,dinvt,dinvp
    double precision,allocatable :: update_value(:,:,:,:)   ! parameter update value 更新参数的变化量

    ! model update (optimization)
    double precision :: stepsize    ! stepsize
    integer :: stIter, edIter  ! max iteration step
    double precision :: old_obj, obj, obj_region, obj_region_l, obj_tele, obj_tele_l, tmp_obj ! 优化目标函数 objective function
    double precision,allocatable :: weight_st(:),weight_ev(:),weight_ev_tele(:) !weight function
    double precision :: dis_ij, dis_zero     !weight function
    logical :: weight_sou_bool=.false. , weight_rec_bool= .false.

    ! temp var
    double precision :: depth,lat,lon
    double precision :: tmp_dp,tmp_int

    ! input and output file
    character(Len=80) :: filename,form,form2,TTfield_path
    character(Len=80) :: input_path,output_path
    logical :: isfile
    logical :: isread_model
    character(Len=80) :: read_model_name

    ! loop index
    integer :: iir,iit,iip,idi,idj,idk,i,j,k,iter,igrid

    ! computing time 
    double precision :: time_begin,time_end,alpha
    
    ! mpi parameter
    integer :: ierr,myid,nproc,tag,istat(mpi_status_size)
    integer,allocatable :: my_tele_source_range(:,:),my_region_source_range(:,:) ! 负载平衡，当前 processor 负责的 最小最大震源序号
    integer :: iproc  ! proc 循环index
    double precision,allocatable :: dp_3d(:,:,:), dp_2d(:,:), dp_1d(:,:), dp_0d


! ############################# input parameters ##############################
    ! this can be included in a input_parameter file later

    ! 正演网格 forward mesh
    nr = 101; nt = 9; np = 101     ! number of grid for forward modeling
    rr1 = 5550; rr2 = 6400          ! 6371 - 5571; 800 km depth
    tt1 = (0.0-1.0)/180*pi; tt2 = (0.0+1.0)/180*pi
    pp1 = (110.0-1.5)/180*pi; pp2 = (130.0+1.5)/180*pi

    ! 反演网格
    ninvr = 15; ninvt = 3; ninvp = 19   ! inversion grid size (coarse)
    invr1 = rr1; invr2 = rr2
    invt1 = tt1 - 1.0/180*pi; invt2 = tt2 + 1.0/180*pi
    invp1 = pp1 - 1.5/180*pi; invp2 = pp2 + 1.5/180*pi
    ngrid = 5
    nk = 3          ! number of typies of kernel
    allocate(inv_dep(ninvr))
    inv_dep = [-25,0,25,50,75,100,150,200,300,400,500,600,700,800,900]


    ! 迭代进程 iterative control
    ! ############### 迭代中途开始反演，以下每个参数都值得仔细考虑和修改
    ! ############### allow starting from one iteration at any step. Please check all following parameters before start the inversion
    stIter = 1; edIter = 50     ! 迭代次数
    isread_model = .false.
    read_model_name = 'ega5/output/model_step100'
    stepsize = 0.01 ! 初始迭代步长


    ! 反演相关参数 parameters for inversion
    max_deg = 2.5/180.0*pi  ! the distance threshold of stations for tele seismic inversion
    slow_p=0.04; ani_p=0.00  ! 检测板异常大小   amplitude of checkerboard anomaly

    region_ratio = 1.0  ! 远近震比例    ratio for regional misfit
    tele_ratio = 1.0    ! ratio for teleseismic misfit
    is_region = .false. ! invert regional tarveltime data or not
    is_tele = .true.    ! invert teleseismic tarveltime data or not

    input_path = 'ega5/input'   ! input and output path
    output_path = 'ega5/output'
    

! ############################# Preparation  #############################

    ! -------------------- iteration check -------------------
    ! 如果当前目录下有生成的文件，暂停程序 if the velocity model file exists, pause
    select case (stIter)
    case (1:9)
        write(form,'(I1)') stIter
    case (10:99)
        write(form,'(I2)') stIter
    case (100:999)
        write(form,'(I3)') stIter
    end select
    filename=trim(output_path)//'/model_step'//trim(form)
    inquire(file=filename,exist = isfile)
    if (isfile) then ! velocity model file exists, pause
        print *, 'file exists, please check the output folder'
        pause
    end if

    ! -------------------- mpi initialization -------------------
    call mpi_init(ierr)

    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    call CPU_TIME(time_begin)   ! program start

    tag = 1 ! communcation tag
    
! ################### read sources and receivers 读取地震台站数据 ###################
    open(10,file=trim(input_path)//'/STATIONS')
    read(10,*) nst
    allocate(rst(nst),tst(nst),pst(nst),weight_st(nst))
    do i=1,nst
        read(10,*) depth,lat,lon
        rst(i) = Earth_Radius-depth
        tst(i) = lat/180.0*pi 
        pst(i) = lon/180.0*pi 
    end do
    close(10)
    
    if (is_region) then
        open(10,file=trim(input_path)//'/EARTHQUAKES_region')
        read(10,*) nev
        allocate(rev(nev),tev(nev),pev(nev),weight_ev(nev))
        do i=1,nev
            read(10,*) depth,lat,lon
            rev(i) = Earth_Radius-depth
            tev(i) = lat/180.0*pi 
            pev(i) = lon/180.0*pi 
        end do
        close(10)
    end if

    if (is_tele) then
        open(10,file=trim(input_path)//'/EARTHQUAKES_tele')
        read(10,*) nev_tele
        allocate(rev_tele(nev_tele),tev_tele(nev_tele),pev_tele(nev_tele),weight_ev_tele(nev_tele))
        do i=1,nev_tele
            read(10,*) depth,lat,lon
            rev_tele(i) = Earth_Radius-depth
            tev_tele(i) = lat/180.0*pi 
            pev_tele(i) = lon/180.0*pi 
        end do
        close(10)
    end if
    ! ---- 数据记录 data record (earthquakes and stations) ----
    if (myid .eq. 0) then
        open(100,file=trim(output_path)//'/STATIONS')
        do i=1,nst
            write(100,'(f13.7,f13.7,f13.7)') &
            & Earth_Radius-rst(i),tst(i)*180.0/pi,pst(i)*180.0/pi
        end do
        close(100); 
        if (is_region) then
            open(101,file=trim(output_path)//'/EARTHQUAKES_region')
            do i=1,nev
                write(101,'(f13.7,f13.7,f13.7)') &
                & Earth_Radius-rev(i),tev(i)*180.0/pi,pev(i)*180.0/pi
            end do     
            close(101);
        end if
        if (is_tele) then
            open(102,file=trim(output_path)//'/EARTHQUAKES_tele')
            do i=1,nev_tele
                write(102,'(f13.7,f13.7,f13.7)') &
                & Earth_Radius-rev_tele(i),tev_tele(i)*180.0/pi,pev_tele(i)*180.0/pi
            end do
            close(102)
        end if
    end if
    
    
! #########################  生成网格 build mesh grid  #####################################
    ! -------- forward modeling grid ----------
    dr=(rr2-rr1)/(nr-1); dt=(tt2-tt1)/(nt-1); dp=(pp2-pp1)/(np-1)
    allocate(rr(nr),tt(nt),pp(np))
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
    allocate(invr(ngrid,ninvr),invt(ngrid,ninvt),invp(ngrid,ninvp))


    do igrid=1,ngrid
        do iir=1,ninvr
            if (iir .eq. 1) then
                dinvr = inv_dep(ninvr)-inv_dep(ninvr-1)
            else
                dinvr = inv_dep(ninvr+1-iir+1)-inv_dep(ninvr+1-iir)
            end if
            invr(igrid,iir)=Earth_Radius-inv_dep(ninvr+1-iir)-(igrid-1)*dinvr/ngrid
        end do
        do iit=1,ninvt
            invt(igrid,iit)=invt1+(iit-1)*dinvt-(igrid-1)*dinvt/ngrid
        end do
        do iip=1,ninvp
            invp(igrid,iip)=invp1+(iip-1)*dinvp-(igrid-1)*dinvp/ngrid
        end do
    end do
    
    ! ------------ 数据记录 data record (multi grid for inversion) -------------
    if (myid==0) then
        open(100,file=trim(output_path)//'/multigrid')
        do igrid=1,Ngrid
            write(100,'(i4,i4,i4,i4)') &
                & igrid,ninvr,ninvt,ninvp
            do iir=1,ninvr
                do iit=1,ninvt
                    do iip=1,ninvp
                        write(100,'(f13.7,f13.7,f13.7)') &
                            & Earth_Radius-invr(igrid,iir),invt(igrid,iit)*180.0/pi,invp(igrid,iip)*180.0/pi
                    end do
                end do
            end do
        end do
        close(100)
    end if

! #########################  构建背景模型 build the synthetic and true model #####################################
    ! ----------- allocation ---------------
    allocate(xi(nr,nt,np),eta(nr,nt,np),zeta(nr,nt,np))   
    allocate(a(nr,nt,np),b(nr,nt,np),c(nr,nt,np),f(nr,nt,np),fun(nr,nt,np))  
    allocate(xiT(nr,nt,np),etaT(nr,nt,np),zetaT(nr,nt,np))  
    allocate(aT(nr,nt,np),bT(nr,nt,np),cT(nr,nt,np),fT(nr,nt,np),funT(nr,nt,np)) 


    filename = trim(input_path)//'/prem_model.txt'
    call Read_1D_velocity_3d(filename,rr,nr,nt,np,fun)

    do iir=1,nr
        do iit=1,nt
            do iip=1,np  

                ! ! ega5 test velocity
                ! if (rr(iir)>6351) then      ! 6371 - 6351   20 km
                !     fun(iir,iit,iip) = 1.0/(5.8+(6371-rr(iir))/20.0*0.7)
                ! elseif (rr(iir)>6336) then  ! 6351 - 6336   15 km
                !     fun(iir,iit,iip) = 1.0/(6.5+(6351-rr(iir))/15.0*0.6)
                ! elseif (rr(iir)>5961) then  ! 6351 - 6336   15 km
                !     fun(iir,iit,iip) = 1.0/(8.0+(6336-rr(iir))/375.0*1)
                ! else
                !     fun(iir,iit,iip) = 1.0/9.0
                ! end if


                ! synthetic (initial) model
                eta(iir,iit,iip)=0.0;
                xi(iir,iit,iip)=0.0;
                zeta(iir,iit,iip)=gamma*sqrt(eta(iir,iit,iip)**2+xi(iir,iit,iip)**2)
                a(iir,iit,iip)=1.0+2*zeta(iir,iit,iip);
                b(iir,iit,iip)=1.0-2*xi(iir,iit,iip);
                c(iir,iit,iip)=1.0+2*xi(iir,iit,iip);
                f(iir,iit,iip)=-2*eta(iir,iit,iip);

                ! true (target) model
                if (pp(iip)>=110.0/180*pi .and. pp(iip)<=130.0/180*pi .and. &
                  & Earth_Radius-rr(iir)>=0.0      .and. Earth_Radius-rr(iir)<=100.0 ) then
                    sigma1 = sin(7*pi*(pp(iip)-110.0/180*pi)/(20.0/180*pi))* &
                           & sin(1*pi*(rr(iir)-6271)/50.0)
                elseif (pp(iip)>=110.0/180*pi .and. pp(iip)<=130.0/180*pi .and. &
                  & Earth_Radius-rr(iir)>=100.0      .and. Earth_Radius-rr(iir)<=200.0 ) then
                    sigma1 = sin(7*pi*(pp(iip)-110.0/180*pi)/(20.0/180*pi))* &
                           & sin(pi+1*pi*(rr(iir)-6171)/100.0)
                elseif (pp(iip)>=110.0/180*pi .and. pp(iip)<=130.0/180*pi .and. &
                  & Earth_Radius-rr(iir)>=200.0      .and. Earth_Radius-rr(iir)<=800.0 ) then
                    sigma1 = sin(7*pi*(pp(iip)-110.0/180*pi)/(20.0/180*pi))* &
                            & sin(1*pi*(rr(iir)-5571)/200.0)           
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
            end do
        end do
    end do
    
    ! 从外部读取模型文件，一般用于从迭代中途开始
    ! if we start the inversion from an input model:
    if (isread_model) then
        call read_model_3d(read_model_name,nr,nt,np,fun,xi,eta,zeta,a,b,c,f)
    else
        ! ---- 数据记录 data record (true and initial model)----
        if (myid .eq. 0) then
            open(100,file=trim(output_path)//'/model_true')
            open(101,file=trim(output_path)//'/model_init')
                do iir=1,nr
                    do iit=1,nt
                        do iip=1,np
                            write(100,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)')  &
                                & Earth_Radius-rr(iir),tt(iit)*180.0/pi,pp(iip)*180.0/pi, &
                                & funT(iir,iit,iip),xiT(iir,iit,iip),etaT(iir,iit,iip)
                            write(101,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)')  &
                                & Earth_Radius-rr(iir),tt(iit)*180.0/pi,pp(iip)*180.0/pi, &
                                & fun(iir,iit,iip),xi(iir,iit,iip),eta(iir,iit,iip)
                        end do
                    end do
                end do
            close(100);close(101)
        end if
    end if

! ######################  震源并行, 负载平衡 load balancing ##################################
    if (is_region) then
        allocate(my_region_source_range(nproc,3))
        ! 分配，每个线程计算那些震源对应的到时
        do iproc=1,nproc
            ! 前 k 个 processor 负责的总震源数，应该接近 震源总数的第k个nproc等分点
            my_region_source_range(iproc,1) = (nst*(iproc-1)/nproc)+1
            my_region_source_range(iproc,2) = (nst*(iproc)/nproc)
            my_region_source_range(iproc,3) = my_region_source_range(iproc,2)-my_region_source_range(iproc,1)+1
        end do
    end if
    if (is_tele) then
        allocate(my_tele_source_range(nproc,3))
        ! 分配，每个线程计算那些震源对应的到时
        do iproc=1,nproc
            ! 前 k 个 processor 负责的总震源数，应该接近 震源总数的第k个nproc等分点
            my_tele_source_range(iproc,1) = (nev_tele*(iproc-1)/nproc)+1
            my_tele_source_range(iproc,2) = (nev_tele*(iproc)/nproc)
            my_tele_source_range(iproc,3) = my_tele_source_range(iproc,2)-my_tele_source_range(iproc,1)+1
        end do
    end if
! ######################  计算远震波到达区域的边界 calculate the traveltime of teleseismic events on the boundary #####################

    if (is_tele) then
        allocate(bd_N(nev_tele,nr,np),bd_S(nev_tele,nr,np),bd_W(nev_tele,nr,nt), &
            & bd_E(nev_tele,nr,nt),bd_B(nev_tele,nr,np),bd_T(nev_tele,nr,np))
        allocate(bd_N_l(nev_tele,nr,np),bd_S_l(nev_tele,nr,np),bd_W_l(nev_tele,nr,nt), &
            & bd_E_l(nev_tele,nr,nt),bd_B_l(nev_tele,nr,np),bd_T_l(nev_tele,nr,np))
        allocate(isbd(nev_tele,5))
        bd_N = 0.0; bd_S = 0.0; bd_W = 0.0; bd_E = 0.0; bd_B = 0.0; bd_T = 0.0; isbd = .false.
        bd_N_l = 0.0; bd_S_l = 0.0; bd_W_l = 0.0; bd_E_l = 0.0; bd_B_l = 0.0; bd_T_l = 0.0;
        filename = trim(input_path)//'/prem_model.txt'
        TTfield_path = trim(output_path)//'/tele_traveltime_field'

        ! calculate the boundary on each processor
        do i = my_tele_source_range(myid+1,1),my_tele_source_range(myid+1,2)      
            print *, '----------------- calculating tele 2d timetable ... evid: ',i,'----------------'
            call eikonal_boundary_traveltime(rr,tt,pp,nr,nt,np,rev_tele(i),tev_tele(i),pev_tele(i) & 
                & ,filename,TTfield_path,isbd(i,:), &
                & bd_N_l(i,:,:),bd_S_l(i,:,:),bd_W_l(i,:,:),bd_E_l(i,:,:),bd_B_l(i,:,:),bd_T_l(i,:,:))
        end do

        
        ! processor communication
        ! 数据通信 步骤1, 其他线程传输给主线程
        call mpi_barrier(mpi_comm_world,ierr)
        if (myid .eq. 0) then   ! 负责接受数据 myid .eq. 0, receiver data
            do iproc=1,nproc-1  
                call mpi_recv(isbd(my_tele_source_range(iproc+1,1):my_tele_source_range(iproc+1,2),:), &
                    & my_tele_source_range(iproc+1,3)*5,mpi_logical,iproc,tag-1,mpi_comm_world,istat,ierr)
            end do
        else    ! 负责发送数据    myid .ne. 0, send data
            call mpi_send(isbd(my_tele_source_range(myid+1,1):my_tele_source_range(myid+1,2),:), &    
                & my_tele_source_range(myid+1,3)*5,mpi_logical,0,tag-1,mpi_comm_world,ierr)                
        end if
        call mpi_allreduce(bd_N_l,bd_N,nev_tele*nr*np,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        call mpi_allreduce(bd_S_l,bd_S,nev_tele*nr*np,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        call mpi_allreduce(bd_W_l,bd_W,nev_tele*nr*nt,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        call mpi_allreduce(bd_E_l,bd_E,nev_tele*nr*nt,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        call mpi_allreduce(bd_B_l,bd_B,nev_tele*nt*np,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        call mpi_allreduce(bd_T_l,bd_T,nev_tele*nt*np,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
        
        ! 数据通信 步骤2, 主线程将更新后的level数据广播到其他线程
        call mpi_bcast(isbd,nev_tele*5,mpi_logical,0,mpi_comm_world,ierr) 
    end if


! #########################  计算真实区域地震到时 regional event true traveltime  ######################
    call mpi_barrier(mpi_comm_world,ierr)
    allocate(TableT(nr,nt,np),uT(nr,nt,np)) ! 初始化走时场 initial the true timetable

    if (is_region) then
        allocate(TtimeT(nev,nst),TtimeT_l(nev,nst)) ! 初始化区域地震真实到时表 initiate true traveltimes 
        TtimeT = 0.0; TtimeT_l = 0.0
        do j=my_region_source_range(myid+1,1),my_region_source_range(myid+1,2)  ! loop STATIONS
            ! --------  求解程函方程 solve eikonal equations ------------
            print *, '----------------- calculating regional true timetable ... evid: ',j,'----------------'
    
            call FSM_WENO3_PS_sphe_3d_mul(rr,tt,pp,nr,nt,np,aT,bT,cT,fT,TableT,funT,rst(j),tst(j),pst(j),uT)
            ! call FSM_WENO3_PS_sphe_3d_mul_mpi(rr,tt,pp,nr,nt,np,aT,bT,cT,fT,TableT,funT,rst(j),tst(j),pst(j),uT)
            ! --------  得到真实走时 calculate true arrival time --------
            do i=1,nev
                call Linear_Interp_3D(rr,tt,pp,TableT,nr,nt,np,rev(i),tev(i),pev(i),TtimeT_l(i,j))
            end do
        end do

        ! processor communication  (true traveltime at stations)
        ! 数据通信 归约
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_allreduce(TtimeT_l,TtimeT,nev*nst,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)


        ! data record (traveltime from source to receiver)
        if ((myid.eq.0) .and. .true.) then
            open(100,file=trim(output_path)//'/traveltime_region')
            do i=1,nev
                do j=1,nst               
                    write(100,'(i4,f8.3,f8.3,f8.3,i4,f8.3,f8.3,f8.3,f8.3)') &
                    & i,Earth_Radius-rev(i),tev(i)/pi*180,pev(i)/pi*180, &
                    & j,Earth_Radius-rst(j),tst(j)/pi*180,pst(j)/pi*180, TtimeT(i,j)
                end do
            end do
            close(100)
        end if
    end if
! #########################  计算真实远震到时 tele seismic event true traveltime  ######################
    if (is_tele) then
        allocate(Ttime_teleT(nev_tele,nst),Ttime_teleT_l(nev_tele,nst));  ! 初始化远震真实到时表 initiate true traveltimes 
        Ttime_teleT = 0.0; Ttime_teleT_l = 0.0 

        do i=my_tele_source_range(myid+1,1),my_tele_source_range(myid+1,2)   ! loop ev_tele
            ! --------  求解远震程函方程 solve eikonal equations ------------
            print *, '----------------- calculating tele true timetable ... evid: ',i,'----------------'
            call FSM_WENO3_tele_sphe_3d_ver2(rr,tt,pp,nr,nt,np,aT,bT,cT,fT,funT,isbd(i,:) &
                                        & ,bd_N(i,:,:),bd_S(i,:,:),bd_W(i,:,:),bd_E(i,:,:),bd_B(i,:,:),TableT)
            ! --------  得到真实走时 calculate true arrival time --------
            do j=1,nst
                call Linear_Interp_3D(rr,tt,pp,TableT,nr,nt,np,rst(j),tst(j),pst(j),Ttime_teleT_l(i,j))
            end do
        end do

        ! processor communication  (true traveltime at stations)
        ! 数据通信 归约
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_allreduce(Ttime_teleT_l,Ttime_teleT,nev_tele*nst,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

        ! data record (traveltime from source to receiver)
        if ((myid.eq.0) .and. .true.) then
            open(100,file=trim(output_path)//'/traveltime_tele')
            do i=1,nev_tele
                do j=1,nst
                    write(100,'(i4,f8.3,f8.3,f8.3,i4,f8.3,f8.3,f8.3,f8.3)') &
                    & i,Earth_Radius-rev_tele(i),tev_tele(i)/pi*180,pev_tele(i)/pi*180, &
                    & j,Earth_Radius-rst(j),tst(j)/pi*180,pst(j)/pi*180, Ttime_teleT(i,j) 
                end do
            end do
            close(100)
        end if
    end if
! ########################  反演开始  inversion start ##########################
    ! ----- 初始化反演参数 initiate inversion parameters ---- 
    old_obj = 0; obj = 0; obj_region=0; obj_tele=0; obj_tele_l=0;! 目标函数 objective function
    
    if (is_region) then
        allocate(all_Ks_region(nr,nt,np),all_Kxi_region(nr,nt,np),all_Keta_region(nr,nt,np))  
        allocate(Ttime(nev,nst))
        Ttime = 0.0; all_Ks_region = 0.0; all_Kxi_region = 0.0; all_Keta_region = 0.0
        allocate(all_Ks_region_l(nr,nt,np),all_Kxi_region_l(nr,nt,np),all_Keta_region_l(nr,nt,np))
        allocate(Ttime_l(nev,nst))
        Ttime_l = 0.0; all_Ks_region_l = 0.0; all_Kxi_region_l = 0.0; all_Keta_region_l = 0.0
    end if
    if (is_tele) then
        allocate(all_Ks_tele(nr,nt,np),all_Kxi_tele(nr,nt,np),all_Keta_tele(nr,nt,np))
        allocate(Ttime_tele(nev_tele,nst))
        Ttime_tele = 0.0; all_Ks_tele = 0.0; all_Kxi_tele = 0.0; all_Keta_tele = 0.0
        allocate(all_Ks_tele_l(nr,nt,np),all_Kxi_tele_l(nr,nt,np),all_Keta_tele_l(nr,nt,np))
        allocate(Ttime_tele_l(nev_tele,nst))
        Ttime_tele_l = 0.0; all_Ks_tele_l = 0.0; all_Kxi_tele_l = 0.0; all_Keta_tele_l = 0.0
        
    end if
    
    allocate(Ks(nr,nt,np),Kxi(nr,nt,np),Keta(nr,nt,np))
    allocate(Table(nr,nt,np),u(nr,nt,np),TableADJ(nr,nt,np))
    allocate(all_Kernel(nr,nt,np,3),update_value(nr,nt,np,3))
    Ks = 0.0; Kxi = 0.0; Keta = 0.0; Table = 0.0; u = 0.0; TableADJ = 0.0;
    all_Kernel = 0.0; update_value = 0.0;

    if (myid .eq. 0) then
        print *, ' '
        print *, '----------------- inversion start ... ----------------'
        print *, ' '
    end if

    do iter = stIter,edIter     ! loop iterations
        call mpi_barrier(mpi_comm_world,ierr)
        if (myid .eq. 0) then
            print *, '----------------- iteration ',iter,' starting ... ----------------'
        end if

        ! ----- 初始化参数 initiate parameters ------
        obj = 0.0; obj_region=0.0; obj_region_l=0.0; obj_tele=0.0; obj_tele_l=0;
        all_Kernel = 0.0; update_value = 0.0;
        
! ########################  区域地震敏感核计算  calculate sensitivity kernel for regional events ##########################      
        if (is_region) then  
            call mpi_barrier(mpi_comm_world,ierr)
            ! ----- 初始化参数 initiate parameters ------
            all_Ks_region_l = 0; all_Kxi_region_l = 0; all_Keta_region_l = 0;
            all_Ks_region = 0; all_Kxi_region = 0; all_Keta_region = 0;

            do j=my_region_source_range(myid+1,1),my_region_source_range(myid+1,2)  ! loop STATIONS
                ! ########################  计算合成走时场 calculate synthetic timetable ########################
                print *, '----------------- calculating regional synthetic timetable and kernel ... evid: ',j,'----------------'
               
                call FSM_WENO3_PS_sphe_3d_mul(rr,tt,pp,nr,nt,np,a,b,c,f,Table,fun,rst(j),tst(j),pst(j),u)

                ! --------  得到合成走时 calculate synthetic arrival time at receiver --------
                do i=1,nev
                    call Linear_Interp_3D(rr,tt,pp,Table,nr,nt,np,rev(i),tev(i),pev(i),Ttime_l(i,j))
                end do
                
                ! ########################  计算伴随场 calculate adjoint timetable ########################

                ! --------  构造伴随源 build adjoint source (absolute traveltime difference) -------
                allocate(sourceADJ(nev))
                call Adjoint_Source_Dt(Ttime_l(:,j),TtimeT(:,j),nev,sourceADJ,tmp_obj)  ! absolute traveltime difference

                call FSM_O1_Adj_sphe_3d(rr,tt,pp,nr,nt,np,Table,TableADJ,zeta,xi,eta,nev,rev,tev,pev,sourceADJ)

                deallocate(sourceADJ)
                ! ######################## 计算敏感核 calculate sensitivity kernel ########################

                call Sensitivity_Kernel(rr,tt,pp,nr,nt,np,Table,TableADJ,gamma,xi,eta,fun,Ks,Kxi,Keta)

                ! ------- 抹去伴随源处的大函数值 mask the source -------   
                call Kernel_Mask_new(rr,tt,pp,nr,nt,np,Ks,rst(j),tst(j),pst(j)) 
                call Kernel_Mask_new(rr,tt,pp,nr,nt,np,Kxi,rst(j),tst(j),pst(j)) 
                call Kernel_Mask_new(rr,tt,pp,nr,nt,np,Keta,rst(j),tst(j),pst(j)) 

                ! --------- 敏感核,obj叠加 sum kernels and objective function -------
                all_Ks_region_l = all_Ks_region_l + Ks;
                all_Kxi_region_l = all_Kxi_region_l + Kxi;
                all_Keta_region_l = all_Keta_region_l + Keta;
                obj_region_l = obj_region_l + tmp_obj

                ! ----- data recored (local event kernel) ------
                if (.false. .and. myid .eq. 0) then
                    select case (j)
                    case (1:9)
                        write(form2,'(I1)') j
                    case (10:99)
                        write(form2,'(I2)') j
                    end select
                    filename=trim(output_path)//'/region_event_kernel_ev'//trim(form2)
                    open(100*(myid+1),file=trim(filename))
                    tmp_dp = 0
                    do iir=1,nr
                        do iit=1,nt
                            do iip=1,np
                                tmp_dp=max(tmp_dp,abs(Ks(iir,iit,iip)))
                                tmp_dp=max(tmp_dp,abs(Kxi(iir,iit,iip)))
                                tmp_dp=max(tmp_dp,abs(Keta(iir,iit,iip)))
                            end do
                        end do
                    end do
                    do iir=1,nr
                        do iit=1,nt
                            do iip=1,np
                                write(100,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)') &
                                & Earth_Radius-rr(iir),tt(iit)*180/pi,pp(iip)*180/pi, &
                                & Ks(iir,iit,iip)/tmp_dp,Kxi(iir,iit,iip)/tmp_dp,Keta(iir,iit,iip)/tmp_dp
                            end do
                        end do
                    end do
                    close(100*(myid+1))
                end if
            end do
    
            ! processor communication  (true traveltime at stations)
            ! 数据通信 归约，把所有进程的kernel和obj求和，传输给所有进程
            call mpi_barrier(mpi_comm_world,ierr)
            call mpi_allreduce(all_Ks_region_l,all_Ks_region,nr*nt*np,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
            call mpi_allreduce(all_Kxi_region_l,all_Kxi_region,nr*nt*np,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
            call mpi_allreduce(all_Keta_region_l,all_Keta_region,nr*nt*np,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
            call mpi_allreduce(obj_region_l,obj_region,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
            call mpi_allreduce(Ttime_l,Ttime,nev*nst,mpi_double_precision,mpi_max,mpi_comm_world,ierr)
  


            ! ----- data recored (regional traveltime misfit at each iteration) ------
            if (myid .eq. 0) then
                select case (iter)
                case (1:9)
                    write(form,'(I1)') iter
                case (10:99)
                    write(form,'(I2)') iter
                case (100:999)
                    write(form,'(I3)') iter
                end select
                filename=trim(output_path)//'/misfit_step'//trim(form)//'_region'
                open(998, file=trim(filename))
                do i=1,nev
                    do j=1,nst     
                        write(998,'(f10.5,f10.5,f10.5)') &
                            & Ttime(i,j)-TtimeT(i,j), Ttime(i,j), TtimeT(i,j)
                    end do
                end do
                close(998)
            end if

            ! ----- data recored (region misfit kernel) ------
            if (.false. .and. myid .eq. 0) then
                filename=trim(output_path)//'/region_misfit_kernel'
                open(100,file=trim(filename))
                tmp_dp = 0
                do iir=1,nr
                    do iit=1,nt
                        do iip=1,np
                            tmp_dp=max(tmp_dp,abs(all_Ks_region(iir,iit,iip)))
                            tmp_dp=max(tmp_dp,abs(all_Kxi_region(iir,iit,iip)))
                            tmp_dp=max(tmp_dp,abs(all_Keta_region(iir,iit,iip)))
                        end do
                    end do
                end do
                do iir=1,nr
                    do iit=1,nt
                        do iip=1,np
                            write(100,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)') &
                            & Earth_Radius-rr(iir),tt(iit)*180/pi,pp(iip)*180/pi, &
                            & all_Ks_region(iir,iit,iip)/tmp_dp,all_Kxi_region(iir,iit,iip)/tmp_dp, &
                            & all_Keta_region(iir,iit,iip)/tmp_dp
                        end do
                    end do
                end do
                close(100)
            end if

        end if
! ########################  远震敏感核计算  calculate sensitivity kernel for teleseismic events ##########################  
        if (is_tele) then  
            call mpi_barrier(mpi_comm_world,ierr)
            ! ----- 初始化参数 initiate parameters ------
            all_Ks_tele_l = 0; all_Kxi_tele_l = 0; all_Keta_tele_l = 0;
            all_Ks_tele = 0; all_Kxi_tele = 0; all_Keta_tele = 0;

            do i=my_tele_source_range(myid+1,1),my_tele_source_range(myid+1,2)
                ! ########################  计算合成走时场 calculate synthetic timetable ########################
                ! --------  求解远震程函方程 solve eikonal equations ------------
                print *, '----------------- calculating tele synthetic timetable and kernel ... evid: ',i,' ----------------'
                
                call FSM_WENO3_tele_sphe_3d_ver2(rr,tt,pp,nr,nt,np,a,b,c,f,fun,isbd(i,:) &
                                                & ,bd_N(i,:,:),bd_S(i,:,:),bd_W(i,:,:),bd_E(i,:,:),bd_B(i,:,:),Table)
                ! --------  得到真实走时 calculate true arrival time --------
                do j=1,nst
                    call Linear_Interp_3D(rr,tt,pp,Table,nr,nt,np,rst(j),tst(j),pst(j),Ttime_tele_l(i,j))
                end do

                ! ########################  计算伴随场 calculate adjoint timetable ########################

                allocate(sourceADJ(nst))
                call Adjoint_Source_DDt(Ttime_tele_l(i,:),Ttime_teleT(i,:),tst,pst,nst,max_deg,sourceADJ,tmp_obj)  ! absolute traveltime difference

                call FSM_O1_Adj_tele_sphe_3d(rr,tt,pp,nr,nt,np,Table,TableADJ,zeta,xi,eta,Nst,rst,tst,pst,sourceADJ,isbd)
                deallocate(sourceADJ)

                ! ########################  计算敏感核 calculate sensitivity kernel ########################

                call Sensitivity_Kernel(rr,tt,pp,nr,nt,np,Table,TableADJ,gamma,xi,eta,fun,Ks,Kxi,Keta)
                
                ! --------- 敏感核,obj叠加 sum kernels and objective function -------
                all_Ks_tele_l = all_Ks_tele_l + Ks;
                all_Kxi_tele_l = all_Kxi_tele_l + Kxi;
                all_Keta_tele_l = all_Keta_tele_l + Keta;
                obj_tele_l = obj_tele_l + tmp_obj
                
                ! ----- data recored (tele event kernel) ------
                if (.false. .and. myid .eq. 0) then
                    select case (i)
                    case (1:9)
                        write(form2,'(I1)') i
                    case (10:99)
                        write(form2,'(I2)') i
                    end select
                    filename=trim(output_path)//'/tele_event_kernel_ev'//trim(form2)
                    open(100*(myid+1),file=trim(filename))
                    tmp_dp = 0
                    do iir=1,nr
                        do iit=1,nt
                            do iip=1,np
                                tmp_dp=max(tmp_dp,abs(Ks(iir,iit,iip)))
                                tmp_dp=max(tmp_dp,abs(Kxi(iir,iit,iip)))
                                tmp_dp=max(tmp_dp,abs(Keta(iir,iit,iip)))
                            end do
                        end do
                    end do
                    do iir=1,nr
                        do iit=1,nt
                            do iip=1,np
                                write(100,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)') &
                                & Earth_Radius-rr(iir),tt(iit)*180/pi,pp(iip)*180/pi, &
                                & Ks(iir,iit,iip)/tmp_dp,Kxi(iir,iit,iip)/tmp_dp,Keta(iir,iit,iip)/tmp_dp
                            end do
                        end do
                    end do
                    close(100*(myid+1))
                end if
            end do

            ! processor communication  (true traveltime at stations)
            ! 数据通信 归约，把所有进程的kernel和obj求和，传输给所有进程
            ! sum kernels and objective functions
            call mpi_barrier(mpi_comm_world,ierr)
            call mpi_allreduce(all_Ks_tele_l,all_Ks_tele,nr*nt*np,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
            call mpi_allreduce(all_Kxi_tele_l,all_Kxi_tele,nr*nt*np,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
            call mpi_allreduce(all_Keta_tele_l,all_Keta_tele,nr*nt*np,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
            call mpi_allreduce(obj_tele_l,obj_tele,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
            call mpi_allreduce(Ttime_tele_l,Ttime_tele,nev_tele*nst,mpi_double_precision,mpi_max,mpi_comm_world,ierr)
          



            ! ----- data recored (teleseismic traveltime misfit at each iteration) ------
            if (myid .eq. 0) then
                select case (iter)
                case (1:9)
                    write(form,'(I1)') iter
                case (10:99)
                    write(form,'(I2)') iter
                case (100:999)
                    write(form,'(I3)') iter
                end select
                filename=trim(output_path)//'/misfit_step'//trim(form)//'_tele'
                open(998, file=trim(filename))
                do i=1,nev_tele
                    do j=1,nst-1    
                        do k=j+1,nst
                            call Epicenter_Distance_sphere(tst(j),pst(j),tst(k),pst(k),deg)
                            if (deg < max_deg) then
                                write(998,'(f10.5,f10.5,f10.5,f10.5,f10.5)') &
                                    & (Ttime_tele(i,j)-Ttime_tele(i,k))-(Ttime_teleT(i,j)-Ttime_teleT(i,k)), &
                                    & Ttime_tele(i,j), Ttime_tele(i,k), Ttime_teleT(i,j), Ttime_teleT(i,k)
                            end if
                        end do
                    end do
                end do
                close(998)
            end if

            ! ----- data recored (tele misfit kernel) ------
            if (.false. .and. myid .eq. 0) then
                filename=trim(output_path)//'/tele_misfit_kernel'
                open(100,file=trim(filename))
                tmp_dp = 0
                do iir=1,nr
                    do iit=1,nt
                        do iip=1,np
                            tmp_dp=max(tmp_dp,abs(all_Ks_tele(iir,iit,iip)))
                            tmp_dp=max(tmp_dp,abs(all_Kxi_tele(iir,iit,iip)))
                            tmp_dp=max(tmp_dp,abs(all_Keta_tele(iir,iit,iip)))
                        end do
                    end do
                end do
                do iir=1,nr
                    do iit=1,nt
                        do iip=1,np
                            write(100,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)') &
                            & Earth_Radius-rr(iir),tt(iit)*180/pi,pp(iip)*180/pi, &
                            & all_Ks_tele(iir,iit,iip)/tmp_dp,all_Kxi_tele(iir,iit,iip)/tmp_dp, &
                            & all_Keta_tele(iir,iit,iip)/tmp_dp
                        end do
                    end do
                end do
                close(100)
            end if

        end if

! ################################### 模型更新 model update ##########################################
        call mpi_barrier(mpi_comm_world,ierr)
        if (myid .eq. 0) then
            print *, '----------------- model update ... iter: ',iter,' ----------------'
        end if

        ! ---- 更新目标函数 update objective function -------
        obj = region_ratio * obj_region + tele_ratio * obj_tele

        ! ----- data recored (misfit and step size at each iteration) ------
        if ((myid .eq. 0)) then
            if (iter .eq. 1) then
                open(999,file=trim(output_path)//'/obj')
            else
                open(999,file=trim(output_path)//'/obj',access='append')
            end if
            write(999,'(f11.2,f11.2,f11.2,f9.6)') obj, obj_region,obj_tele,stepsize
            close(999)
        end if

        ! --------- 整合kernel (sum kernels) ---------------------
        if (is_region) then
            all_Kernel(:,:,:,1) = all_Kernel(:,:,:,1) + region_ratio * all_Ks_region 
            all_Kernel(:,:,:,2) = all_Kernel(:,:,:,2) + region_ratio * all_Kxi_region 
            all_Kernel(:,:,:,3) = all_Kernel(:,:,:,3) + region_ratio * all_Keta_region 
        end if
        if (is_tele) then
            all_Kernel(:,:,:,1) = all_Kernel(:,:,:,1) + tele_ratio * all_Ks_tele
            all_Kernel(:,:,:,2) = all_Kernel(:,:,:,2) + tele_ratio * all_Kxi_tele
            all_Kernel(:,:,:,3) = all_Kernel(:,:,:,3) + tele_ratio * all_Keta_tele
        end if
        

        ! ------------ multiple grid parameterization ------------
        ! call Parameter_Update_Multigrid(rr,tt,pp,nr,nt,np,all_Kernel,nk, &
        !                                 & invr,invt,invp,ninvr,ninvt,ninvp,ngrid,stepsize,update_value)
        call Parameter_Update_Multigrid_ver2(rr,tt,pp,nr,nt,np,all_Kernel,nk, &
                                        & invr,invt,invp,ninvr,ninvt,ninvp,ngrid,stepsize,update_value)
        do iir=1,nr
            do iit=1,nt
                do iip=1,np
                    ! model parameter update
                    fun(iir,iit,iip) = fun(iir,iit,iip)*(1-update_value(iir,iit,iip,1))
                    xi(iir,iit,iip) = xi(iir,iit,iip)-update_value(iir,iit,iip,2)
                    eta(iir,iit,iip) = eta(iir,iit,iip)-update_value(iir,iit,iip,3)

                    b(iir,iit,iip)=1.0-2*xi(iir,iit,iip);
                    c(iir,iit,iip)=1.0+2*xi(iir,iit,iip);
                    f(iir,iit,iip)=-2*eta(iir,iit,iip);
                end do
            end do
        end do

! ################################## 调整下降步长 modify stepsize ##############################

        if (iter == 1 ) then
            if (myid .eq. 0) then
                write(*,'(a,f9.2)') 'iter 1, obj is', obj
                write(*,'(a,f9.6)') 'iter 1, stepsize is', stepsize
            end if
        elseif (iter >= 2 .and. obj < old_obj) then
            if (myid .eq. 0) then
                write(*,'(a,f9.2,a,f9.2)') 'objective function decreases, from', old_obj, ' to', obj
                write(*,'(a,f9.6)') 'new stepsize is ', stepsize
            end if
        elseif (iter >= 2 .and. obj >= old_obj) then
            stepsize = max(0.003,stepsize*0.97)   ! ega5 会注释掉步长变化
            if (myid .eq. 0) then
                write(*,'(a,f9.2,a,f9.2)') 'objective function increases, from', old_obj, ' to', obj
                write(*,'(a,f9.6)') 'new stepsize is ', stepsize
            end if
        end if

! ################################## 迭代结束 iteration over ##############################
        if (myid .eq. 0) then
            print *, ' '
            print *, '----------------- iteration ',iter,' over ----------------'
            print *, ' '
        end if

        old_obj = obj

        ! ---- data recored (parameters at each iteration) ----
        if ((myid .eq. 0) .and. .true.) then
            select case (iter)
            case (1:9)
                write(form,'(I1)') iter
            case (10:99)
                write(form,'(I2)') iter
            case (100:999)
                write(form,'(I3)') iter
            end select

            filename=trim(output_path)//'/model_step'//trim(form)
            open(100,file=trim(filename))
            do iir=1,nr
                do iit=1,nt
                    do iip=1,np
                        write(100,'(f13.7,f13.7,f13.7,f13.7,f13.7,f13.7)') &
                            & Earth_Radius-rr(iir),tt(iit)*180.0/pi,pp(iip)*180.0/pi, &
                            & fun(iir,iit,iip),xi(iir,iit,iip),eta(iir,iit,iip)
                    end do
                end do
            end do
            close(100)
        end if

        
        if (myid .eq. 0) then
            call CPU_TIME(time_end)
            write(*,'(a,f10.2,a)') 'runtime is: ',time_end - time_begin, 'second'
        end if



    end do

! 储存释放 release space
    if (is_tele) then
        deallocate(rev_tele,tev_tele,pev_tele,weight_ev_tele)
        deallocate(bd_N,bd_S,bd_W,bd_E,bd_B,bd_T,isbd)
        deallocate(Ttime_teleT,Ttime_tele)
        deallocate(Ttime_teleT_l,Ttime_tele_l)
        deallocate(all_Ks_tele,all_Kxi_tele,all_Keta_tele)
        deallocate(all_Ks_tele_l,all_Kxi_tele_l,all_Keta_tele_l)
        deallocate(my_tele_source_range)
    end if
    if (is_region) then
        deallocate(rev,tev,pev,weight_ev)
        deallocate(TtimeT,Ttime)
        deallocate(TtimeT_l,Ttime_l)
        deallocate(all_Ks_region,all_Kxi_region,all_Keta_region)
        deallocate(all_Ks_region_l,all_Kxi_region_l,all_Keta_region_l)
        deallocate(my_region_source_range)
    end if
    
    deallocate(rst,tst,pst,weight_st)
    deallocate(rr,tt,pp,invr,invt,invp)
    deallocate(xi,eta,zeta)   
    deallocate(a,b,c,f,fun)  
    deallocate(xiT,etaT,zetaT)  
    deallocate(aT,bT,cT,fT,funT) 
    deallocate(TableT,uT,Table,u,TableADJ)

    deallocate(Ks,Kxi,Keta)
    deallocate(all_kernel,update_value)

    deallocate(inv_dep)


    call mpi_finalize(ierr)

end program eikonal_2d

