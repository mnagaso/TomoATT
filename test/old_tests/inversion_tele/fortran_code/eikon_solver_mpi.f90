subroutine FSM_WENO3_tele_sphe_3d(rr,tt,pp,nr,nt,np,spha,sphb,sphc,sphf,T,fun,u,ischange)
    ! a -d -e
    ! -d b -f
    ! -e -f c
    integer nr,nt,np
    double precision :: dr,dt,dp,rr(nr),tt(nt),pp(np),a(nr,nt,np),b(nr,nt,np),c(nr,nt,np),f(nr,nt,np)
    double precision :: spha(nr,nt,np),sphb(nr,nt,np),sphc(nr,nt,np),sphf(nr,nt,np)
    double precision :: fun(nr,nt,np),T(nr,nt,np),ischange(nr,nt,np),u(nr,nt,np),xi(nr,nt,np),eta(nr,nt,np)
    double precision :: T_old(nr,nt,np)
    double precision :: sigr,sigt,sigp,coe,pr1,pr2,wr1,wr2,pt1,pt2,wt1,wt2,pp1,pp2,wp1,wp2,tpT
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err
    double precision,parameter :: tol = (10.0)**(-5),eps=10.0**(-12)
    integer,parameter :: MaxIter=500
    integer iter,rdirec,tdirec,pdirec,iir,iit,iip
    double precision :: tmp  
    
    ! ------------------------ 构造网格 ------------------------ 
    dr=rr(2)-rr(1); dt=tt(2)-tt(1); dp=pp(2)-pp(1)

     ! ------------------------ 构造矩阵 build eikonal matrix -------------------------  
    !    a -d -e
    !   -d  b -f
    !   -e -f  c
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                a(iir,iit,iip) = spha(iir,iit,iip)
                b(iir,iit,iip) = sphb(iir,iit,iip)/(rr(iir)**2)
                c(iir,iit,iip) = sphc(iir,iit,iip)/(rr(iir)**2*cos(tt(iit))**2)
                f(iir,iit,iip) = sphf(iir,iit,iip)/(rr(iir)**2*cos(tt(iit)))
            end do
        end do
    end do

    ! 正式迭代，更新 iteration start, update T
    do iter =1,MaxIter
        T_old = T
        L1_dif=0; Linf_dif=0;L1_err=0;Linf_err=0;
        do rdirec = -1,1,2
            do tdirec = -1,1,2
                do pdirec = -1,1,2

                    !x: nr <-> 1, y: nt <-> 1, z: np <-> 1

                    do iir=nint(0.5+nr/2.0+(nr/2.0-0.5)*rdirec),nint(0.5+nr/2.0+(-nr/2.0+0.5)*rdirec),-rdirec
                        do iit=nint(0.5+nt/2.0+(nt/2.0-0.5)*tdirec),nint(0.5+nt/2.0+(-nt/2.0+0.5)*tdirec),-tdirec
                            do iip=nint(0.5+np/2.0+(np/2.0-0.5)*pdirec),nint(0.5+np/2.0+(-np/2.0+0.5)*pdirec),-pdirec
                                

                                if(ischange(iir,iit,iip)==1) then
                                    if( iir>1 .and. iir<nr .and. iit>1 .and. iit<Nt .and. iip>1 .and. iip<Np) then

                                        sigr=sqrt(a(iir,iit,iip)); sigt=sqrt(b(iir,iit,iip)); sigp=sqrt(c(iir,iit,iip))
                                        coe=1.0/((sigr/dr)+(sigt/dt)+(sigp/dp))

                                        ! forward and backward partial derivatives

                                        
                                            pr1=(T(iir,iit,iip)-T(iir-1,iit,iip))/dr;
                                            pr2=(T(iir+1,iit,iip)-T(iir,iit,iip))/dr;                                       
                        
                                            pt1=(T(iir,iit,iip)-T(iir,iit-1,iip))/dt; 
                                            pt2=(T(iir,iit+1,iip)-T(iir,iit,iip))/dt;  
                                    
                                            pp1=(T(iir,iit,iip)-T(iir,iit,iip-1))/dp;                    
                                            pp2=(T(iir,iit,iip+1)-T(iir,iit,iip))/dp; 


                                        ! calculate  LF Hamiltonian 

                                        HT=sqrt( a(iir,iit,iip)*((pr1+pr2)/2)**2 + b(iir,iit,iip)*((pt1+pt2)/2)**2 &
                                        & + c(iir,iit,iip)*((pp1+pp2)/2)**2 -2*f(iir,iit,iip)*(pt1+pt2)/2*(pp1+pp2)/2 )
                                        
                                        ! 更新 update timetable
                                        tpT=coe*(fun(iir,iit,iip)-HT)  &
                                        & +coe*(sigr*(pr2-pr1)/2+sigt*(pt2-pt1)/2+sigp*(pp2-pp1)/2)+T(iir,iit,iip);

                                        !write(*,*) fun(iir,iit,iip),HT

                                    
                                        T(iir,iit,iip) = tpT
                                        

                                    elseif (iir==1) then
                                        T(1,iit,iip) = max(2*T(2,iit,iip)-T(3,iit,iip),T(3,iit,iip))
                                    elseif (iir==nr) then
                                        T(nr,iit,iip) = max(2*T(nr-1,iit,iip)-T(nr-2,iit,iip),T(nr-2,iit,iip))
                                    elseif (iit==1) then        
                                        T(iir,1,iip) = max(2*T(iir,2,iip)-T(iir,3,iip),T(iir,3,iip))
                                    elseif (iit==nt) then        
                                        T(iir,nt,iip) = max(2*T(iir,nt-1,iip)-T(iir,nt-2,iip),T(iir,nt-2,iip))
                                    elseif (iip==1) then        
                                        T(iir,iit,1) = max(2*T(iir,iit,2)-T(iir,iit,3),T(iir,iit,3))
                                    elseif (iip==Np) then 
                                        T(iir,iit,np) = max(2*T(iir,iit,np-1)-T(iir,iit,np-2),T(iir,iit,np-2))

                                    end if
                                end if

                                
                            end do
                        end do
                    end do    

                end do
            end do
        end do

        
        ! 统计误差，判断迭代终止条件
        do iir=1,nr
            do iit=1,nt
                do iip=1,np
                    L1_dif=L1_dif+abs(T(iir,iit,iip)-T_old(iir,iit,iip))
                    Linf_dif=max(Linf_dif,abs(T(iir,iit,iip)-T_old(iir,iit,iip)))
                end do               
            end do
        end do
        
        do iir=1,nr
            do iit=1,nt
                do iip=1,np
                    L1_err=L1_err+abs(u(iir,iit,iip)-T(iir,iit,iip))
                    Linf_err=max(Linf_err,abs(T(iir,iit,iip)-u(iir,iit,iip)))   
                end do            
            end do
        end do
        L1_err=L1_err/((nr)*(nt)*(np))
        L1_dif=L1_dif/((nr)*(nt)*(np))


        ! ################  iteration information  #################
        if (abs(L1_dif)<tol) then
            write(*,*) 'iter ',iter,', T is steadt'
            exit
        else
            write(*,*) 'iter ',iter,', T is changing, L1 dif = ', L1_dif,'L inf dif = ', Linf_dif
        end if

        if (iter==MaxIter) then    
            write(*,*) 'iter ',iter,', max iteration steps'
        end if

        ! ################  solver accuracy  #################
        !write(*,'(a,f15.7,a,f15.7)') 'L_1(T-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err


    end do

end subroutine

! ----------------------- for ega1 ----------------------
subroutine read_boundary(rr,tt,pp,nr,nt,np,T,ischange,isbd,u,surface)
    integer :: nr,nt,np
    double precision :: T(nr,nt,np),rr(nr),tt(nt),pp(np),ischange(nr,nt,np)
    double precision :: u(nr,nt,np),surface(nt,np)
    double precision,parameter :: Radius = 6371.0
    double precision,parameter :: pi=3.14159265358979323846264338327950288
    double precision,parameter :: infT = 2000.0
    logical :: isbd(6),isN,isS,isW,isE,isT,isB  ! N S W E T B

    !  6 个方向的边界
    double precision,allocatable :: bd_W(:,:), bd_E(:,:), bd_N(:,:), bd_S(:,:), bd_T(:,:), bd_B(:,:)
    double precision,allocatable :: bd_rr(:), bd_tt(:), bd_pp(:)
    double precision :: bd_W_interp(nr,nt), bd_E_interp(nr,nt), bd_N_interp(nr,np), bd_S_interp(nr,np)
    double precision :: bd_T_interp(nt,np), bd_B_interp(nt,np)

    ! 边界的index
    integer :: ni,nj
    double precision  :: x1,x2,y1,y2
    integer :: i,j
    double precision  :: read_tmp

    isN=isbd(1);isS=isbd(2)
    isW=isbd(3);isE=isbd(4)
    isT=isbd(5);isB=isbd(6)
    

    ! #---------------------------   西 W 边界  ---------------------------#  
    open(10, file='ega1/input/boundary_condition_W')
    read(10, *) ni,nj,x1,x2,y1,y2
    allocate(bd_W(ni,nj),bd_rr(ni),bd_tt(nj))
    do i=1,ni
        do j=1,nj
            read(10, *) read_tmp, read_tmp, read_tmp, bd_W(i,j), read_tmp, read_tmp
        end do
    end do
    x1 = Radius - x1
    x2 = Radius - x2
    y1 = y1/180*pi
    y2 = y2/180*pi
    
    bd_dr = (x2-x1)/(ni-1)
    bd_dt = (y2-y1)/(nj-1)
    do i=1,ni
        bd_rr(i) = x1 + (i-1)*bd_dr
    end do
    do j=1,nj
        bd_tt(j) = y1 + (j-1)*bd_dt
    end do

    ! 远震到时 插值到 计算区域边界
    call boundary_interp(bd_rr,bd_tt,ni,nj,bd_W,rr,tt,nr,nt,bd_W_interp)
    
    ! 赋值给 走时场T
    do i=1,nr
        do j=1,nt
            if (isW) then
                if (bd_W_interp(i,j)>0) then
                    T(i,j,1) = bd_W_interp(i,j)
                    u(i,j,1) = bd_W_interp(i,j)
                    ischange(i,j,1) = 0
                else
                    T(i,j,1) = infT
                    ischange(i,j,1) = 1
                end if
            else
                if (bd_W_interp(i,j)>0) then
                    u(i,j,1) = bd_W_interp(i,j)
                end if
            end if

        end do
    end do

    ! 释放动态数组
    deallocate(bd_rr,bd_tt,bd_W)
    

    
    ! #---------------------------   东 E 边界  ---------------------------#  
    open(10, file='ega1/input/boundary_condition_E')
    read(10, *) ni,nj,x1,x2,y1,y2
    allocate(bd_E(ni,nj),bd_rr(ni),bd_tt(nj))
    do i=1,ni
        do j=1,nj
            read(10, *) read_tmp, read_tmp, read_tmp, bd_E(i,j), read_tmp, read_tmp
        end do
    end do
    x1 = Radius - x1
    x2 = Radius - x2
    y1 = y1/180*pi
    y2 = y2/180*pi
    
    bd_dr = (x2-x1)/(ni-1)
    bd_dt = (y2-y1)/(nj-1)
    do i=1,ni
        bd_rr(i) = x1 + (i-1)*bd_dr
    end do
    do j=1,nj
        bd_tt(j) = y1 + (j-1)*bd_dt
    end do

    ! 远震到时 插值到 计算区域边界
    call boundary_interp(bd_rr,bd_tt,ni,nj,bd_E,rr,tt,nr,nt,bd_E_interp)
    
    ! 赋值给 走时场T
    do i=1,nr
        do j=1,nt
            if (isE) then
                if (bd_E_interp(i,j)>0) then
                    T(i,j,np) = bd_E_interp(i,j)
                    u(i,j,np) = bd_E_interp(i,j)
                    ischange(i,j,np) = 0
                else
                    T(i,j,np) = infT
                    ischange(i,j,np) = 1
                end if
            else
                if (bd_E_interp(i,j)>0) then
                    u(i,j,np) = bd_E_interp(i,j)
                end if
            end if
        end do
    end do

    ! 释放动态数组
    deallocate(bd_rr,bd_tt,bd_E)
    

    
    ! #---------------------------   北 N 边界  ---------------------------#  
    open(10, file='ega1/input/boundary_condition_N')
    read(10, *) ni,nj,x1,x2,y1,y2
    allocate(bd_N(ni,nj),bd_rr(ni),bd_pp(nj))
    do i=1,ni
        do j=1,nj
            read(10, *) read_tmp, read_tmp, read_tmp, bd_N(i,j), read_tmp, read_tmp
        end do
    end do
    x1 = Radius - x1
    x2 = Radius - x2
    y1 = y1/180*pi
    y2 = y2/180*pi
    
    bd_dr = (x2-x1)/(ni-1)
    bd_dp = (y2-y1)/(nj-1)
    do i=1,ni
        bd_rr(i) = x1 + (i-1)*bd_dr
    end do
    do j=1,nj
        bd_pp(j) = y1 + (j-1)*bd_dp
    end do

    ! 远震到时 插值到 计算区域边界
    call boundary_interp(bd_rr,bd_pp,ni,nj,bd_N,rr,pp,nr,np,bd_N_interp)
    
    ! 赋值给 走时场T
    do i=1,nr
        do j=1,np
            if (isN) then
                if (bd_N_interp(i,j)>0) then
                    T(i,nt,j) = bd_N_interp(i,j)
                    u(i,nt,j) = bd_N_interp(i,j)
                    ischange(i,nt,j) = 0
                else
                    T(i,nt,j) = infT
                    ischange(i,nt,j) = 1
                end if
            else
                if (bd_N_interp(i,j)>0) then
                    u(i,nt,j) = bd_N_interp(i,j)
                end if
            end if
        end do
    end do

    ! 释放动态数组
    deallocate(bd_rr,bd_pp,bd_N)
    

    
    ! #---------------------------   南 S 边界  ---------------------------#  
    open(10, file='ega1/input/boundary_condition_S')
    read(10, *) ni,nj,x1,x2,y1,y2
    allocate(bd_S(ni,nj),bd_rr(ni),bd_pp(nj))
    do i=1,ni
        do j=1,nj
            read(10, *) read_tmp, read_tmp, read_tmp, bd_S(i,j), read_tmp, read_tmp
        end do
    end do
    x1 = Radius - x1
    x2 = Radius - x2
    y1 = y1/180*pi
    y2 = y2/180*pi
    
    bd_dr = (x2-x1)/(ni-1)
    bd_dp = (y2-y1)/(nj-1)
    do i=1,ni
        bd_rr(i) = x1 + (i-1)*bd_dr
    end do
    do j=1,nj
        bd_pp(j) = y1 + (j-1)*bd_dp
    end do

    ! 远震到时 插值到 计算区域边界
    call boundary_interp(bd_rr,bd_pp,ni,nj,bd_S,rr,pp,nr,np,bd_S_interp)
    
    ! 赋值给 走时场T
    do i=1,nr
        do j=1,np
            if (isS) then
                if (bd_S_interp(i,j)>0) then
                    T(i,1,j) = bd_S_interp(i,j)
                    u(i,1,j) = bd_S_interp(i,j)
                    ischange(i,1,j) = 0
                else
                    T(i,1,j) = infT
                    ischange(i,1,j) = 1
                end if
            else
                if (bd_S_interp(i,j)>0) then
                    u(i,1,j) = bd_S_interp(i,j)
                end if
            end if
        end do
    end do

    ! 释放动态数组
    deallocate(bd_rr,bd_pp,bd_S)
    

    ! #---------------------------   顶 Top 边界  ---------------------------#  
    open(10, file='ega1/input/boundary_condition_T')
    read(10, *) ni,nj,x1,x2,y1,y2
    allocate(bd_T(ni,nj),bd_tt(ni),bd_pp(nj))
    do i=1,ni
        do j=1,nj
            read(10, *) read_tmp, read_tmp, read_tmp, bd_T(i,j), read_tmp, read_tmp
        end do
    end do
    x1 = x1/180*pi
    x2 = x2/180*pi
    y1 = y1/180*pi
    y2 = y2/180*pi
    
    bd_dt = (x2-x1)/(ni-1)
    bd_dp = (y2-y1)/(nj-1)
    do i=1,ni
        bd_tt(i) = x1 + (i-1)*bd_dt
    end do
    do j=1,nj
        bd_pp(j) = y1 + (j-1)*bd_dp
    end do

    ! 远震到时 插值到 计算区域边界
    call boundary_interp(bd_tt,bd_pp,ni,nj,bd_T,tt,pp,nt,np,bd_T_interp)
    
    ! 赋值给 走时场T
    do i=1,nt
        do j=1,np
            if (isT) then
                if (bd_T_interp(i,j)>0) then
                    T(nr,i,j) = bd_T_interp(i,j)
                    ischange(nr,i,j) = 0
                    surface(i,j) = bd_T_interp(i,j)
                else
                    T(nr,i,j) = infT
                    ischange(nr,i,j) = 1
                end if
            else
                if (bd_T_interp(i,j)>0) then
                    surface(i,j) = bd_T_interp(i,j)
                end if
            end if
        end do
    end do

    ! 释放动态数组
    deallocate(bd_tt,bd_pp,bd_T)

    
    ! #---------------------------   底 Bottom 边界  ---------------------------#  
    open(10, file='ega1/input/boundary_condition_B')
    read(10, *) ni,nj,x1,x2,y1,y2
    allocate(bd_B(ni,nj),bd_tt(ni),bd_pp(nj))
    do i=1,ni
        do j=1,nj
            read(10, *) read_tmp, read_tmp, read_tmp, bd_B(i,j), read_tmp, read_tmp
        end do
    end do
    x1 = x1/180*pi
    x2 = x2/180*pi
    y1 = y1/180*pi
    y2 = y2/180*pi
    
    bd_dt = (x2-x1)/(ni-1)
    bd_dp = (y2-y1)/(nj-1)
    do i=1,ni
        bd_tt(i) = x1 + (i-1)*bd_dt
    end do
    do j=1,nj
        bd_pp(j) = y1 + (j-1)*bd_dp
    end do

    ! 远震到时 插值到 计算区域边界
    call boundary_interp(bd_tt,bd_pp,ni,nj,bd_B,tt,pp,nt,np,bd_B_interp)
    
    ! 赋值给 走时场T
    do i=1,nt
        do j=1,np
            if (isB) then
                if (bd_T_interp(i,j)>0) then
                    T(1,i,j) = bd_B_interp(i,j)
                    u(1,i,j) = bd_B_interp(i,j)
                    ischange(1,i,j) = 0
                else
                    T(1,i,j) = infT
                    ischange(1,i,j) = 1
                end if
            else
                if (bd_T_interp(i,j)>0) then
                    u(1,i,j) = bd_B_interp(i,j)
                end if
            end if
        end do
    end do

    ! 释放动态数组
    deallocate(bd_tt,bd_pp,bd_B)
    
end subroutine

subroutine boundary_interp(xx1,yy1,nx1,ny1,val1,xx2,yy2,nx2,ny2,val2)
    ! 将 1 套网格的值 插值到 2套网格上，允许外插
    integer :: nx1,ny1,nx2,ny2
    double precision :: xx1(nx1),yy1(ny1),xx2(nx2),yy2(ny2)
    double precision :: dx1,dy1,val1(nx1,ny1),val2(nx2,ny2)
    double precision :: r1,r2
    integer i,j,idx0,idy0
    dx1 = xx1(2)-xx1(1)
    dy1 = yy1(2)-yy1(1)



    do i=1,nx2
        do j=1,ny2
            idx0 = floor((xx2(i)-xx1(1))/dx1)+1
            idy0 = floor((yy2(j)-yy1(1))/dy1)+1

            if (idx0<1 .or. idx0>nx1-1 .or. idy0<1 .or. idy0>ny1-1) then
            ! 允许外插
                idx0 = min(nx1-1,max(1,idx0))
                idy0 = min(ny1-1,max(1,idy0)) 
            end if
            
            
            r1 = min(1.0, (xx2(i)-xx1(idx0))/dx1 )
            r2 = min(1.0, (yy2(j)-yy1(idy0))/dy1 )
            
            val2(i,j) = (1-r1)*(1-r2)*val1(idx0,idy0) + (1-r1)*(r2)*val1(idx0,idy0+1) &
                    & + (r1)*(1-r2)*val1(idx0+1,idy0) + (r1)*(r2)*val1(idx0+1,idy0+1)
                  
        end do
    end do
end subroutine

subroutine tele_error_estimation(rr,tt,pp,nr,nt,np,T,filename,dep,error)
    integer :: nr,nt,np
    double precision :: T(nr,nt,np),rr(nr),tt(nt),pp(np)
    double precision :: error(nt,np,3),dep
    double precision,parameter :: Radius = 6371.0
    double precision,parameter :: pi=3.14159265358979323846264338327950288
    character(Len=80) :: filename
    
    double precision,allocatable :: tauP_solution(:,:),tauP_tt(:),tauP_pp(:)
    double precision :: tauP_dt,tauP_dp,tauP_solution_interp(nt,np)

    double precision :: L1_err, Linf_err

    ! 边界的index
    integer :: ni,nj
    double precision  :: x1,x2,y1,y2,r1
    integer :: i,j,idx0,iit,iip
    double precision  :: read_tmp
    
    dr = rr(2)-rr(1)


    ! 读取指定文件
    open(10, file=filename)
    read(10, *) ni,nj,x1,x2,y1,y2
    allocate(tauP_solution(ni,nj),tauP_tt(ni),tauP_pp(nj))
    

    do i=1,ni
        do j=1,nj
            read(10, *) read_tmp, read_tmp, read_tmp, tauP_solution(i,j), read_tmp, read_tmp
        end do
    end do
    x1 = x1/180*pi
    x2 = x2/180*pi
    y1 = y1/180*pi
    y2 = y2/180*pi
    
    tauP_dt = (x2-x1)/(ni-1)
    tauP_dp = (y2-y1)/(nj-1)
    do i=1,ni
        tauP_tt(i) = x1 + (i-1)*tauP_dt
    end do
    do j=1,nj
        tauP_pp(j) = y1 + (j-1)*tauP_dp
    end do

    ! 远震到时 插值到 计算区域
    call boundary_interp(tauP_tt,tauP_pp,ni,nj,tauP_solution,tt,pp,nt,np,tauP_solution_interp)
    
    ! 走时场 T 插值到指定深度
    idx0 = floor(((Radius-dep) -rr(1))/dr)+1
    r1 = min(1.0,((Radius-dep)-rr(idx0))/dr )
    L1_err=0; Linf_err=0

    do iit=1,nt
        do iip=1,np
            error(iit,iip,1) = tauP_solution_interp(iit,iip)
            error(iit,iip,2) = (1-r1)*T(idx0,iit,iip)+r1*T(idx0+1,iit,iip)       
            error(iit,iip,3) = error(iit,iip,1) - error(iit,iip,2)
            L1_err = L1_err + abs(error(iit,iip,3))
            Linf_err = max(Linf_err, abs(error(iit,iip,3)))
        end do
    end do
    L1_err = L1_err/(nt*np)
    ! 释放动态数组
    deallocate(tauP_solution,tauP_tt,tauP_pp)
    
    write(*,'(a,f5.1,a,es10.3,a,es10.3,a)') 'L1 and Linf errors at depth ', dep, 'km are ', L1_err, ' s and ', Linf_err,' s.'


end subroutine
! ----------------------- end for ega1 ----------------------


! ----------------------- for ega2 ----------------------
subroutine load_boundary(rr,tt,pp,nr,nt,np,T,ischange,taup_fn,isbd,tele_location)
    integer :: nr,nt,np
    double precision :: T(nr,nt,np),rr(nr),tt(nt),pp(np),ischange(nr,nt,np)
    double precision,parameter :: Radius = 6371.0,pi=3.14159265358979323846264338327950288
    logical :: isbd(6),isN,isS,isW,isE,isT,isB  ! N S W E T B
    double precision :: tele_location(3)
    
    ! tauP database
    character(Len=80) :: taup_fn
    integer :: taup_nsrc_dep,taup_nrec_dep,taup_ndegree
    double precision,allocatable :: taup_time(:,:,:),taup_src_dep(:),taup_rec_dep(:),taup_degree(:)
    double precision :: tmp_read    ! useless for read

    ! boundary points
    double precision :: epicenter_dis(nt,np),lat1,lon1,lat2,lon2
    double precision :: r1,r2,r3    ! linear interpolation 
    double precision :: src,rec,dis
    integer :: id1,id2,id3
    ! otheres 
    integer :: i,j,k    !iterative index


    ! 读取 tauP 数据库 三列分别是 震源深度，台站深度，和震中距
    open(10, file=taup_fn)
    read(10, *) taup_nsrc_dep,taup_nrec_dep,taup_ndegree

    allocate(taup_time(taup_nsrc_dep,taup_nrec_dep,taup_ndegree))
    allocate(taup_src_dep(taup_nsrc_dep),taup_rec_dep(taup_nrec_dep),taup_degree(taup_ndegree))
     
    read(10, *) taup_src_dep
    read(10, *) taup_rec_dep
    read(10, *) taup_degree

    do i=1,taup_nsrc_dep
        do j=1,taup_nrec_dep
            do k=1,taup_ndegree
                read(10, *) tmp_read,tmp_read,tmp_read,taup_time(i,j,k)
            end do
        end do
    end do
    close(10)
    ! 插值边界

    isN=isbd(1);isS=isbd(2)
    isW=isbd(3);isE=isbd(4)
    isT=isbd(5);isB=isbd(6)

    ! 统一计算震中距
    do i=1,nt
        do j=1,np
            lat1 = tt(i); lon1 = pp(j)
            lat2 = tele_location(2); lon2 = tele_location(3); 
            epicenter_dis(i,j)=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))/pi*180
        end do
    end do


    ! 北边界 North 
    if (isN) then
        do i=1,nr
            do j=1,np
                ! 震源位置: tele_location(1), 台站位置 Radius-rr(i)，震中距：epicenter_dis(nt,j),
                src = tele_location(1)
                rec = Radius-rr(i)
                dis = epicenter_dis(nt,j)

                call taup_data_interp(taup_src_dep,taup_rec_dep,taup_degree, &
                            & taup_time,taup_nsrc_dep,taup_nrec_dep,taup_ndegree,src,rec,dis,T(i,nt,j))  
                ischange(i,nt,j) = 0           
            end do
        end do
    end if

    ! 南边界 South
    if (isS) then
        do i=1,nr
            do j=1,np
                ! 震源位置: tele_location(1), 台站位置 Radius-rr(i)，震中距：epicenter_dis(1,j),
                src = tele_location(1)
                rec = Radius-rr(i)
                dis = epicenter_dis(1,j)

                call taup_data_interp(taup_src_dep,taup_rec_dep,taup_degree, &
                            & taup_time,taup_nsrc_dep,taup_nrec_dep,taup_ndegree,src,rec,dis,T(i,1,j))             
                ischange(i,1,j) = 0    
            end do
        end do
    end if

    ! 西边界 West
    if (isW) then
        do i=1,nr
            do j=1,nt
                ! 震源位置: tele_location(1), 台站位置 Radius-rr(i)，震中距：epicenter_dis(j,1),
                src = tele_location(1)
                rec = Radius-rr(i)
                dis = epicenter_dis(j,1)

                call taup_data_interp(taup_src_dep,taup_rec_dep,taup_degree, &
                            & taup_time,taup_nsrc_dep,taup_nrec_dep,taup_ndegree,src,rec,dis,T(i,j,1))    
                ischange(i,j,1) = 0          
            end do
        end do
    end if

    ! 东边界 East
    if (isE) then
        do i=1,nr
            do j=1,nt
                ! 震源位置: tele_location(1), 台站位置 Radius-rr(i)，震中距：epicenter_dis(j,np),
                src = tele_location(1)
                rec = Radius-rr(i)
                dis = epicenter_dis(j,np)

                call taup_data_interp(taup_src_dep,taup_rec_dep,taup_degree, &
                            & taup_time,taup_nsrc_dep,taup_nrec_dep,taup_ndegree,src,rec,dis,T(i,j,np))     
                ischange(i,j,np) = 0         
            end do
        end do
    end if

    ! 上边界 top (一般是不插值的，遂略去)

    ! 下边界 bottom
    if (isB) then
        do i=1,nt
            do j=1,np
                ! 震源位置: tele_location(1), 台站位置 Radius-rr(1)，震中距：epicenter_dis(i,j),
                src = tele_location(1)
                rec = Radius-rr(1)
                dis = epicenter_dis(i,j)

                call taup_data_interp(taup_src_dep,taup_rec_dep,taup_degree, &
                            & taup_time,taup_nsrc_dep,taup_nrec_dep,taup_ndegree,src,rec,dis,T(1,i,j))  
                ischange(1,i,j) = 0            
            end do
        end do
    end if


    deallocate(taup_time,taup_src_dep,taup_rec_dep,taup_degree)

end subroutine

subroutine taup_data_interp(src,rec,deg,time,ns,nr,nd,s,r,d,v)
    integer :: ns,nr,nd
    double precision :: time(ns,nr,nd),src(ns),rec(nr),deg(nd)
    double precision :: s,r,d,v

    integer :: id1,id2,id3,k
    double precision ::  r1,r2,r3
    ! 震源插值
    ! do k=1,ns-1
    !     if (src(k)<= s .and. src(k+1) > s ) then
    !         id1=k
    !         exit
    !     elseif (src(1) > s) then
    !         id1=1
    !         exit
    !     elseif (src(ns) <= s) then
    !         id1=ns-1
    !         exit
    !     end if
    ! end do
    ! r1 = (s - src(id1)) / (src(id1+1)-src(id1))
    
    id1 = 1

    ! 台站插值 (允许外插)
    do k=1,nr-1
        if (rec(k)<= r .and. rec(k+1) > r ) then
            id2=k
            exit
        elseif (rec(1) > r) then
            id2=1
            exit
        elseif (rec(nr) <= r) then
            id2=nr-1
            exit
        end if
    end do
    r2 = (r - rec(id2)) / (rec(id2+1)-rec(id2))

    ! 震中距插值 (允许外插)
    do k=1,nd-1
        if (deg(k)<= d .and. deg(k+1) > d ) then
            id3=k
            exit
        elseif (deg(1) > d) then
            id3=1
            exit
        elseif (deg(nd) <= d) then
            id3=nd-1
            exit
        end if
    end do
    r3 = (d - deg(id3)) / (deg(id3+1)-deg(id3))

    ! v = (1-r1)*(1-r2)*(1-r3)*time(id1,id2,id3) + (1-r1)*(1-r2)*(r3)*time(id1,id2,id3+1) &
    ! & + (1-r1)*(r2)*(1-r3)*time(id1,id2+1,id3) + (1-r1)*(r2)*(r3)*time(id1,id2+1,id3+1) &
    ! & + (r1)*(1-r2)*(1-r3)*time(id1+1,id2,id3) + (r1)*(1-r2)*(r3)*time(id1+1,id2,id3+1) &
    ! & + (r1)*(r2)*(1-r3)*time(id1+1,id2+1,id3) + (r1)*(r2)*(r3)*time(id1+1,id2+1,id3+1)

    v = (1-r2)*(1-r3)*time(id1,id2,id3) + (1-r2)*(r3)*time(id1,id2,id3+1) &
    & + (r2)*(1-r3)*time(id1,id2+1,id3) + (r2)*(r3)*time(id1,id2+1,id3+1)
    

end subroutine

subroutine tele_error_estimation_ver2(rr,tt,pp,nr,nt,np,T,taup_fn,tele_location,dep,error)
    integer :: nr,nt,np
    double precision :: T(nr,nt,np),rr(nr),tt(nt),pp(np)
    double precision :: error(nt,np,3),dep,tele_location(3)
    double precision,parameter :: Radius = 6371.0
    double precision,parameter :: pi=3.14159265358979323846264338327950288
    
    ! tauP database
    character(Len=80) :: taup_fn
    integer :: taup_nsrc_dep,taup_nrec_dep,taup_ndegree
    double precision,allocatable :: taup_time(:,:,:),taup_src_dep(:),taup_rec_dep(:),taup_degree(:)
    double precision :: tmp_read    ! useless for read

    double precision :: epicenter_dis(nt,np),lat1,lon1,lat2,lon2
    double precision :: src,rec,dis,traveltime

    double precision :: L1_err, Linf_err



    ! 边界的index
    integer :: ni,nj
    double precision  :: x1,x2,y1,y2,r1
    integer :: i,j,k,idx0,iit,iip
    double precision  :: read_tmp
    
    dr = rr(2)-rr(1)


    ! 读取 tauP 数据库 三列分别是 震源深度，台站深度，和震中距
    open(10, file=taup_fn)
    read(10, *) taup_nsrc_dep,taup_nrec_dep,taup_ndegree

    allocate(taup_time(taup_nsrc_dep,taup_nrec_dep,taup_ndegree))
    allocate(taup_src_dep(taup_nsrc_dep),taup_rec_dep(taup_nrec_dep),taup_degree(taup_ndegree))
     
    read(10, *) taup_src_dep
    read(10, *) taup_rec_dep
    read(10, *) taup_degree
    
    do i=1,taup_nsrc_dep
        do j=1,taup_nrec_dep
            do k=1,taup_ndegree
                read(10, *) tmp_read,tmp_read,tmp_read,taup_time(i,j,k)
            end do
        end do
    end do
    close(10)
    ! 统一计算震中距
    do i=1,nt
        do j=1,np
            lat1 = tt(i); lon1 = pp(j)
            lat2 = tele_location(2); lon2 = tele_location(3); 
            epicenter_dis(i,j)=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))/pi*180
        end do
    end do

    ! 指定深度 tauP 走时
    idx0 = floor(((Radius-dep) -rr(1))/dr)+1
    r1 = min(1.0,((Radius-dep)-rr(idx0))/dr )
    L1_err=0; Linf_err=0

    do iit=1,nt
        do iip=1,np
            ! 震源位置: tele_location(1), 台站位置 Radius-dep，震中距：epicenter_dis(iit,iip),
            src = tele_location(1)
            rec = dep
            dis = epicenter_dis(iit,iip)

            call taup_data_interp(taup_src_dep,taup_rec_dep,taup_degree, &
                        & taup_time,taup_nsrc_dep,taup_nrec_dep,taup_ndegree,src,rec,dis,traveltime)      

            error(iit,iip,1) = traveltime
            error(iit,iip,2) = (1-r1)*T(idx0,iit,iip)+r1*T(idx0+1,iit,iip)       
            error(iit,iip,3) = error(iit,iip,1) - error(iit,iip,2)
            L1_err = L1_err + abs(error(iit,iip,3))
            Linf_err = max(Linf_err, abs(error(iit,iip,3)))

        end do
    end do
    L1_err = L1_err/(nt*np)

    ! 释放动态数组
    deallocate(taup_time,taup_src_dep,taup_rec_dep,taup_degree)
    
    write(*,'(a,f5.1,a,es10.3,a,es10.3,a)') 'L1 and Linf errors at depth ', dep, 'km are ', L1_err, ' s and ', Linf_err,' s.'


end subroutine
! ----------------------- end for ega2 ----------------------


! ----------------------- for ega3 ----------------------

subroutine FSM_WENO3_PS_sphe_2d(rr,tt,nr,nt,spha,sphb,T,fun,r0,t0,u)
    integer :: nr,nt
    double precision :: spha(nr,nt),sphb(nr,nt)
    double precision :: rr(nr),tt(nt),a(nr,nt),b(nr,nt),T(nr,nt),fun(nr,nt),r0,t0,ischange(nr,nt),u(nr,nt)
    double precision :: dr,dt,T0v(nr,nt),T0r(nr,nt),T0t(nr,nt),a0,b0,fun0
    double precision :: tau(nr,nt),px1,px2,py1,py2,tpT,Htau,wx1,wx2,wy1,wy2
    integer :: iir,iit,iter,xdirec,ydirec
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err,tau_old(nr,nt),sigr,sigt
    integer,parameter :: MaxIter = 2000
    double precision,parameter :: tol=(10.0)**(-4),eps=(10.0)**(-12)

    ! ------------------------ 构造网格 ------------------------ 
    dr=rr(2)-rr(1); dt=tt(2)-tt(1);

    ! ------------------------ 构造 T0 ------------------------ 

    ! 震源处参数离散化
    idr0=floor((r0-rr(1))/dr+1); idt0=floor((t0-tt(1))/dt+1); 
    r1 = min(1.0,(r0-rr(idr0))/dr); r2 = min(1.0,(t0-tt(idt0))/dt); 

    do iir=1,nr
        do iit=1,nt
            a(iir,iit) = spha(iir,iit)
            b(iir,iit) = sphb(iir,iit)/(rr(iir)**2)
        end do
    end do

    a0=(1-r1)*(1-r2)*a(idr0,idt0)+(1-r1)*r2*a(idr0,idt0+1) &
    & +r1*(1-r2)*a(idr0+1,idt0)+r1*r2*a(idr0+1,idt0+1) 

    b0=(1-r1)*(1-r2)*b(idr0,idt0)+(1-r1)*r2*b(idr0,idt0+1) &
    & +r1*(1-r2)*b(idr0+1,idt0)+r1*r2*b(idr0+1,idt0+1) 

    fun0=(1-r1)*(1-r2)*fun(idr0,idt0)+(1-r1)*r2*fun(idr0,idt0+1) &
    & +r1*(1-r2)*fun(idr0+1,idt0)+r1*r2*fun(idr0+1,idt0+1) 

    do iir=1,nr
        do iit=1,nt
            T0v(iir,iit) = fun0*sqrt((1.0/a0)*(rr(iir)-r0)**2 + 1.0/b0*(tt(iit)-t0)**2)
            if (T0v(iir,iit) .eq. 0) then
                T0r(iir,iit) = 0
                T0t(iir,iit) = 0
            else
                T0r(iir,iit) = fun0**2*(1.0/a0*(rr(iir)-r0))/T0v(iir,iit)
                T0t(iir,iit) = fun0**2*(1.0/b0*(tt(iit)-y0))/T0v(iir,iit)
            end if

            if ( abs((rr(iir)-r0)/dr)<=2 .and. abs((tt(iit)-t0)/dt)<=2) then
                tau(iir,iit) = 1  !震源周围几个点，直接认为是常速度结构，给出解析解，即和T0相等
                ischange(iir,iit)=0
                if (iir==1 .or. iir==nr .or. iit==1 .or. iit==nt) then
                    write(*,*) 'source on the boundary, mesh error'
                    pause
                end if
            else
                tau(iir,iit) = 1
                ischange(iir,iit)=1
            end if

        end do
    end do

    ! step 2, solve Tau, H(tau) = a tau_x^2+ b tau_y^2 + (2aTx-2cTy) tau_x + (2bTy-2cTx) tau_y
    !                           -2c tau_x tau_y + (aTx^2+bTy^2-2cTxTy) = f^2

    do iter =1,MaxIter
        L1_dif=10000; Linf_dif=10000;L1_err=0;Linf_err=0;
        tau_old = tau
        do xdirec = -1,1,2
            do ydirec = -1,1,2
                ! iter 1 x: 1 -> nr, y: 1 -> nt
                ! iter 2 x: 1 -> nr, y: nt -> 1 
                ! iter 3 x: nr -> 1, y: 1 -> nt 
                ! iter 4 x: nr -> 1, y: nt -> 1 
                do iir=nint(0.5+nr/2.0+(nr/2.0-1.5)*xdirec),nint(0.5+nr/2.0+(-nr/2.0+1.5)*xdirec),-xdirec
                    do iit=nint(0.5+nt/2.0+(nt/2.0-1.5)*ydirec),nint(0.5+nt/2.0+(-nt/2.0+1.5)*ydirec),-ydirec
                        if(ischange(iir,iit)==1) then
                            
                            sigr=sqrt(a(iir,iit))*T0v(iir,iit);sigt=sqrt(b(iir,iit))*T0v(iir,iit)
                            coe=1.0/((sigr/dr)+(sigt/dt)) 

                            if(iir==2) then
                                px1=(tau(iir,iit)-tau(iir-1,iit))/dr; 
                                wx2=1.0/(1+2*((eps+(tau(iir,iit)-2*tau(iir+1,iit)+tau(iir+2,iit))**2)/ &
                                            & (eps+(tau(iir-1,iit)-2*tau(iir,iit)+tau(iir+1,iit))**2))**2)
                                px2=(1-wx2)*(tau(iir+1,iit)-tau(iir-1,iit))/2/dr+ &
                                      & wx2*(-3*tau(iir,iit)+4*tau(iir+1,iit)-tau(iir+2,iit))/2/dr;
                            elseif (iir==nr-1) then
                                wx1=1.0/(1+2*((eps+(tau(iir,iit)-2*tau(iir-1,iit)+tau(iir-2,iit))**2)/ &
                                            & (eps+(tau(iir+1,iit)-2*tau(iir,iit)+tau(iir-1,iit))**2))**2)
                                px1=(1-wx1)*(tau(iir+1,iit)-tau(iir-1,iit))/2/dr+ &
                                      & wx1*(3*tau(iir,iit)-4*tau(iir-1,iit)+tau(iir-2,iit))/2/dr;
                                px2=(tau(iir+1,iit)-tau(iir,iit))/dr; 
                            else
                                wx1=1.0/(1.0+2*((eps+(tau(iir,iit)-2*tau(iir-1,iit)+tau(iir-2,iit))**2)/ &
                                            & (eps+(tau(iir+1,iit)-2*tau(iir,iit)+tau(iir-1,iit))**2))**2)
                                px1=(1.0-wx1)*(tau(iir+1,iit)-tau(iir-1,iit))/2/dr+ &
                                      & wx1*(3*tau(iir,iit)-4*tau(iir-1,iit)+tau(iir-2,iit))/2/dr;
                                wx2=1.0/(1.0+2*((eps+(tau(iir,iit)-2*tau(iir+1,iit)+tau(iir+2,iit))**2)/ &
                                            & (eps+(tau(iir-1,iit)-2*tau(iir,iit)+tau(iir+1,iit))**2))**2)
                                px2=(1.0-wx2)*(tau(iir+1,iit)-tau(iir-1,iit))/2/dr+ &
                                      & wx2*(-3*tau(iir,iit)+4*tau(iir+1,iit)-tau(iir+2,iit))/2/dr;
                            end if


                            if(iit==2) then
                                py1=(tau(iir,iit)-tau(iir,iit-1))/dt; 
                                wy2=1.0/(1+2*((eps+(tau(iir,iit)-2*tau(iir,iit+1)+tau(iir,iit+2))**2)/ &
                                            & (eps+(tau(iir,iit-1)-2*tau(iir,iit)+tau(iir,iit+1))**2))**2)
                                py2=(1-wy2)*(tau(iir,iit+1)-tau(iir,iit-1))/2/dt+ &
                                    & wy2*(-3*tau(iir,iit)+4*tau(iir,iit+1)-tau(iir,iit+2))/2/dt;
                            elseif (iit==nt-1) then
                                wy1=1.0/(1+2*((eps+(tau(iir,iit)-2*tau(iir,iit-1)+tau(iir,iit-2))**2)/ &
                                            & (eps+(tau(iir,iit+1)-2*tau(iir,iit)+tau(iir,iit-1))**2))**2)
                                py1=(1-wy1)*(tau(iir,iit+1)-tau(iir,iit-1))/2/dt+ &
                                    & wy1*(3*tau(iir,iit)-4*tau(iir,iit-1)+tau(iir,iit-2))/2/dt;
                                py2=(tau(iir,iit+1)-tau(iir,iit))/dt;
                            else
                                wy1=1.0/(1+2*((eps+(tau(iir,iit)-2*tau(iir,iit-1)+tau(iir,iit-2))**2)/ &
                                            & (eps+(tau(iir,iit+1)-2*tau(iir,iit)+tau(iir,iit-1))**2))**2)
                                py1=(1-wy1)*(tau(iir,iit+1)-tau(iir,iit-1))/2/dt+ &
                                    & wy1*(3*tau(iir,iit)-4*tau(iir,iit-1)+tau(iir,iit-2))/2/dt;
                                wy2=1.0/(1+2*((eps+(tau(iir,iit)-2*tau(iir,iit+1)+tau(iir,iit+2))**2)/ &
                                            & (eps+(tau(iir,iit-1)-2*tau(iir,iit)+tau(iir,iit+1))**2))**2)
                                py2=(1-wy2)*(tau(iir,iit+1)-tau(iir,iit-1))/2/dt+ &
                                    & wy2*(-3*tau(iir,iit)+4*tau(iir,iit+1)-tau(iir,iit+2))/2/dt;                            
                            end if
                


                            Htau = sqrt( a(iir,iit)*(T0r(iir,iit)*tau(iir,iit)+T0v(iir,iit)*(px1+px2)/2)**2 &
                                     & + b(iir,iit)*(T0t(iir,iit)*tau(iir,iit)+T0v(iir,iit)*(py1+py2)/2)**2 )

                            tpT=coe*(fun(iir,iit)-Htau) + coe*(sigr*(px2-px1)/2+sigt*(py2-py1)/2)+tau(iir,iit);
                                                        
                            tau(iir,iit)=tpT


                        end if

                    end do
                end do
                
                ! 处理边界
                do iit=1,nt
                    tau(1,iit) = max(2*tau(2,iit)-tau(3,iit),tau(3,iit))
                    tau(nr,iit) = max(2*tau(nr-1,iit)-tau(nr-2,iit),tau(nr-2,iit))
                end do
                do iir=1,nr
                    tau(iir,1) = max(2*tau(iir,2)-tau(iir,3),tau(iir,3))
                    tau(iir,nt) = max(2*tau(iir,nt-1)-tau(iir,nt-2),tau(iir,nt-2))
                end do
                
            end do
        end do

        L1_dif=0; Linf_dif=0
        do iir=1,nr
            do iit=1,nt
                L1_dif=L1_dif+abs(tau(iir,iit)-tau_old(iir,iit))
                Linf_dif=max(Linf_dif,abs(tau(iir,iit)-tau_old(iir,iit)))               
            end do
        end do
        L1_dif = L1_dif/(nr*nt)

        do iir=3,nr-2
            do iit=3,nt-2
                L1_err=L1_err+abs(tau(iir,iit)*T0v(iir,iit)-u(iir,iit))
                Linf_err=max(Linf_err,abs(tau(iir,iit)*T0v(iir,iit)-u(iir,iit)))            
            end do
        end do
        L1_err = L1_err/(nr-4)/(nt-4)


        if (abs(L1_dif)<tol .and. Linf_dif<tol) then
            write(*,*) 'iter ',iter,', T is steady'
            exit
        else
        !    write(*,*) 'iter ',iter,', T is changing, continue ... '
        end if

        if (iter==MaxIter) then    
            write(*,*) 'iter ',iter,', max iteration steps'
            exit
        end if
        !write(*,'(a,f10.7,a,f10.7)') 'L_1(Tnew-Told)=',L1_dif,'  L_inf(Tnew-Told)',Linf_dif
        ! write(*,'(a,f15.4,a,f15.4)') 'L_1(T_LF-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err
    end do

    T=tau*T0v

end subroutine

subroutine Read_1D_velocity_2d(filename,rr,nr,nt,fun)
    integer :: nr,nt
    double precision :: fun(nr,nt),rr(nr)
    character(len=80) filename

    ! input velocity file
    integer :: nlayer
    double precision,allocatable :: layer_dep(:), layer_vel(:)
    double precision,parameter :: Radius = 6371.0, small_eps=0.001
    
    ! others
    integer i,iir,iit   ! iterative index
    double precision :: read_tmp,dep,vel    ! read temp variable
    


    open(10,file=trim(filename))
    read(10,*) nlayer
    allocate(layer_dep(nlayer),layer_vel(nlayer))

    do i=1,nlayer
        read(10,*) layer_dep(i),layer_vel(i)
    end do

    do iir=1,nr
        ! 计算 dep 处的速度
        dep = Radius - rr(iir)
        if (dep < layer_dep(1)) then
            vel = layer_vel(1)
        elseif (dep >= layer_dep(nlayer)) then
            vel = layer_vel(nlayer)
        else
            do i=1,nlayer-1
                if ( layer_dep(i)+small_eps <= dep .and. layer_dep(i+1) > dep-small_eps ) then
                    vel = (dep-layer_dep(i))/(layer_dep(i+1)-layer_dep(i))*(layer_vel(i+1)-layer_vel(i))+layer_vel(i)
                    exit
                end if
            end do
        end if

        do iit=1,nt
            fun(iir,iit) = 1.0/vel
        end do

    end do

    close(10)
    deallocate(layer_dep,layer_vel)


end subroutine

subroutine Read_1D_velocity_3d(filename,rr,nr,nt,np,fun)
    integer :: nr,nt,np
    double precision :: fun(nr,nt,np),rr(nr)
    character(len=80) filename

    ! input velocity file
    integer :: nlayer
    double precision,allocatable :: layer_dep(:), layer_vel(:)
    double precision,parameter :: Radius = 6371.0, small_eps=0.001
    
    ! others
    integer i,iir,iit,iip   ! iterative index
    double precision :: read_tmp,dep,vel    ! read temp variable
    


    open(10,file=trim(filename))
    read(10,*) nlayer
    allocate(layer_dep(nlayer),layer_vel(nlayer))

    do i=1,nlayer
        read(10,*) layer_dep(i),layer_vel(i)
    end do

    do iir=1,nr
        ! 计算 dep 处的速度
        dep = Radius - rr(iir)
        if (dep < layer_dep(1)) then
            vel = layer_vel(1)
        elseif (dep >= layer_dep(nlayer)) then
            vel = layer_vel(nlayer)
        else
            do i=1,nlayer-1
                if ( layer_dep(i)+small_eps <= dep .and. layer_dep(i+1) > dep-small_eps ) then
                    vel = (dep-layer_dep(i))/(layer_dep(i+1)-layer_dep(i))*(layer_vel(i+1)-layer_vel(i))+layer_vel(i)
                    exit
                end if
            end do
        end if

        do iit=1,nt
            do iip=1,np
                fun(iir,iit,iip) = 1.0/vel
            end do
        end do

    end do

    close(10)
    deallocate(layer_dep,layer_vel)


end subroutine

subroutine Ray_tracing_sphe_2d(rr,tt,nr,nt,T,recr,rect,srcr,srct,time,filename)
    integer :: nr,nt
    double precision :: rr(nr),tt(nt),T(nr,nt)
    double precision :: srcr,srct,recr,rect,time
    character(len=80) :: filename

    integer :: Ntime,traceid
    double precision,allocatable :: trace(:,:)
    double precision :: dt,grad_r,grad_t,Tgrad_r(nr,nt),Tgrad_t(nr,nt)
    double precision :: tol

    integer :: iir,iit,it   ! 循环 index
    double precision,parameter :: pi=3.14159265358979323846264338327950288

    

    dtime=0.1;
    Ntime = int(time*1.5/dtime)   ! 总时间
    allocate(trace(Ntime,2))

    ! 初始化
    trace(1,1) = recr; trace(1,2) = rect; ! 起始点

    dr = rr(2)-rr(1); dt=tt(2)-tt(1)
    tol = sqrt(dr**2+(dt*srcr)**2)
    

    ! 梯度场
    do iir=2,nr-1
        do iit=2,nt-1
            Tgrad_r(iir,iit) = (T(iir+1,iit) - T(iir-1,iit))/dr/2
            Tgrad_t(iir,iit) = (T(iir,iit+1) - T(iir,iit-1))/dt/2
        end do
    end do
    do iir=1,nr
        Tgrad_t(iir,1) = (T(iir,2) - T(iir,1))/dt
        Tgrad_t(iir,nt) = (T(iir,nt) - T(iir,nt-1))/dt
    end do
    do iit=1,nt
        Tgrad_r(1,iit) = (T(2,iit) - T(1,iit))/dr
        Tgrad_r(nr,iit) = (T(nr,iit) - T(nr-1,iit))/dr
    end do

    traceid=-1

    

    ! 主循环
    do it=1,Ntime-1

        ! 该点位置 trace(it,1), trace(it,2)
        if ( trace(it,1)<=rr(1) .or. trace(it,1)>=rr(nr) .or. trace(it,2)<=tt(1) .or. trace(it,2)>=tt(nt)) then
            print *, 'ray out of region, stop'
            exit
        end if
        
        ! 如果当前点离震源非常近，判定收敛
        if ( (trace(it,1)-srcr)**2+(trace(it,2)-srct)**2*trace(it,1)**2 < tol**2 ) then
            print *, 'ray arrive source'
            trace(it+1,1) = srcr
            trace(it+1,2) = srct
            traceid = it+1
            exit
        end if

        
        ! 计算该点的梯度
        call interp2d(rr,tt,nr,nt,Tgrad_r,trace(it,1),trace(it,2),grad_r)
        call interp2d(rr,tt,nr,nt,Tgrad_t,trace(it,1),trace(it,2),grad_t)

        slowness_squared = grad_r**2 + grad_t**2/trace(it,1)**2
        ! 沿负梯度走 dtime 长度  r = r- dtime*Tr/(s^2); theta = theta - dtime*Ttheta/(r^2 s^2)
        trace(it+1,1) = trace(it,1) - dtime*grad_r / slowness_squared
        trace(it+1,2) = trace(it,2) - dtime*grad_t / slowness_squared / ((trace(it,1)+trace(it+1,1))**2/4)

        !print *, trace(it+1,1),trace(it+1,2),xsou,ysou

    end do

    if (traceid<0) then
        ! 没打到震源，存在问题
        print *, filename, 'not arrive'
        open(10,file=trim(filename)//'_NotArrive')
        do it=1,Ntime
            write(10,*) trace(it,1),trace(it,2)*180.0/pi,dtime*(it-1)
        end do
        close(10)
    else
        ! 达到震源，输出射线路径
        open(10,file=trim(filename))
        do it=1,traceid
            write(10,*) trace(it,1),trace(it,2)*180.0/pi,dtime*(it-1)
        end do
        close(10)
    end if

    

    deallocate(trace)
end subroutine

subroutine interp2d(xx,yy,nx,ny,val,x0,y0,v0)
    integer :: nx,ny,idx0,idy0
    double precision :: xx(nx),yy(ny),val(nx,ny),x0,y0,v0,dx,dy,r1,r2
    integer :: iix,iiy

    dx=xx(2)-xx(1); dy=yy(2)-yy(1)


    if ( x0<=xx(1) .or. x0>=xx(nx) .or. y0<=yy(1) .or. y0>=yy(ny)) then
        print *, x0,xx(1),xx(nx),y0,yy(1),yy(ny)
        print *, 'out of mesh'
        pause
    end if

    idx0 = floor((x0-xx(1))/dx)+1
    idy0 = floor((y0-yy(1))/dy)+1

    r1 = min(1.0, (x0-xx(idx0))/dx )
    r2 = min(1.0, (y0-yy(idy0))/dy )

    v0 = (1-r1)*(1-r2)*val(idx0,idy0) + (1-r1)*r2*val(idx0,idy0+1) &
    & + r1*(1-r2)*val(idx0+1,idy0) + r1*r2*val(idx0+1,idy0+1)


end subroutine

subroutine load_boundary_ver2(rr,tt,pp,nr,nt,np,T,ischange,bd_fn,isbd,tele_location)
    integer :: nr,nt,np
    double precision :: T(nr,nt,np),rr(nr),tt(nt),pp(np),ischange(nr,nt,np)
    double precision,parameter :: Radius = 6371.0,pi=3.14159265358979323846264338327950288
    logical :: isbd(6),isN,isS,isW,isE,isT,isB  ! N S W E T B
    double precision :: tele_location(3)
    
    ! boundary database
    character(Len=80) :: bd_fn
    integer :: bd_nr,bd_nt
    double precision :: bd_srcr
    double precision,allocatable :: bd_time(:,:),bd_r(:),bd_t(:)
    double precision :: tmp_read    ! useless for read

    ! boundary points
    double precision :: epicenter_dis(nt,np),lat1,lon1,lat2,lon2
    double precision :: r1,r2,r3    ! linear interpolation 
    double precision :: rec,dis
    integer :: id1,id2,id3
    ! otheres 
    integer :: i,j,k    !iterative index

    ! 读取 2d eikonal 计算的到时数据
    open(10, file=bd_fn)
    read(10, *) bd_srcr
    if (abs(bd_srcr - tele_location(1))>0.1) then
        print *, 'source depth is not consistent with the boundary data'
        pause
    end if
    
    read(10, *) bd_nr,bd_nt
    
    allocate(bd_time(bd_nr,bd_nt),bd_r(bd_nr),bd_t(bd_nt))
     
    read(10, *) bd_r
    read(10, *) bd_t

    do i=1,bd_nr
        do j=1,bd_nt
            read(10, *) bd_time(i,j)
        end do
    end do
    close(10)
    ! 插值边界

    

    isN=isbd(1);isS=isbd(2)
    isW=isbd(3);isE=isbd(4)
    isT=isbd(5);isB=isbd(6)

    ! 统一计算震中距
    do i=1,nt
        do j=1,np
            lat1 = tt(i); lon1 = pp(j)
            lat2 = tele_location(2); lon2 = tele_location(3); 
            epicenter_dis(i,j)=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))/pi*180
        end do
    end do

    
    ! 北边界 North 
    if (isN) then
        do i=1,nr
            do j=1,np
                ! 震源位置: tele_location(1), 台站位置 Radius-rr(i)，震中距：epicenter_dis(nt,j),
                rec = rr(i)
                dis = epicenter_dis(nt,j)
                
                call interp2d(bd_r,bd_t,bd_nr,bd_nt,bd_time,rec,dis,T(i,nt,j))
                ischange(i,nt,j) = 0           
            end do
        end do
    end if

    ! 南边界 South
    if (isS) then
        do i=1,nr
            do j=1,np
                ! 震源位置: tele_location(1), 台站位置 Radius-rr(i)，震中距：epicenter_dis(1,j),
                rec = rr(i)
                dis = epicenter_dis(1,j)
                call interp2d(bd_r,bd_t,bd_nr,bd_nt,bd_time,rec,dis,T(i,1,j))
                ischange(i,1,j) = 0    
            end do
        end do
    end if

    ! 西边界 West
    if (isW) then
        do i=1,nr
            do j=1,nt
                ! 震源位置: tele_location(1), 台站位置 Radius-rr(i)，震中距：epicenter_dis(j,1),
                rec = rr(i)
                dis = epicenter_dis(j,1)

                call interp2d(bd_r,bd_t,bd_nr,bd_nt,bd_time,rec,dis,T(i,j,1))

                ischange(i,j,1) = 0          
            end do
        end do
    end if

    ! 东边界 East
    if (isE) then
        do i=1,nr
            do j=1,nt
                ! 震源位置: tele_location(1), 台站位置 Radius-rr(i)，震中距：epicenter_dis(j,np),
                rec = rr(i)
                dis = epicenter_dis(j,np)

                call interp2d(bd_r,bd_t,bd_nr,bd_nt,bd_time,rec,dis,T(i,j,np))

                ischange(i,j,np) = 0         
            end do
        end do
    end if

    ! 上边界 top (一般是不插值的，遂略去)

    ! 下边界 bottom
    if (isB) then
        do i=1,nt
            do j=1,np
                ! 震源位置: tele_location(1), 台站位置 Radius-rr(1)，震中距：epicenter_dis(i,j),
                rec = rr(1)
                dis = epicenter_dis(i,j)

                call interp2d(bd_r,bd_t,bd_nr,bd_nt,bd_time,rec,dis,T(1,i,j))

                ischange(1,i,j) = 0            
            end do
        end do
    end if


    deallocate(bd_time,bd_r,bd_t)

end subroutine

subroutine tele_error_estimation_ver3(rr,tt,pp,nr,nt,np,T,bd_fn,tele_location,dep,error)
    integer :: nr,nt,np
    double precision :: T(nr,nt,np),rr(nr),tt(nt),pp(np)
    double precision :: error(nt,np,3),dep,tele_location(3)
    double precision,parameter :: Radius = 6371.0
    double precision,parameter :: pi=3.14159265358979323846264338327950288
    
    ! boundary database
    character(Len=80) :: bd_fn
    integer :: bd_nr,bd_nt
    double precision :: bd_srcr
    double precision,allocatable :: bd_time(:,:),bd_r(:),bd_t(:)
    double precision :: tmp_read    ! useless for read


    double precision :: epicenter_dis(nt,np),lat1,lon1,lat2,lon2
    double precision :: src,rec,dis,traveltime
    double precision :: L1_err, Linf_err



    ! 边界的index
    integer :: ni,nj
    double precision  :: x1,x2,y1,y2,r1
    integer :: i,j,k,idx0,iit,iip
    double precision  :: read_tmp
    
    dr = rr(2)-rr(1)


   ! 读取 2d eikonal 计算的到时数据
    open(10, file=bd_fn)
    read(10, *) bd_srcr
    if (abs(bd_srcr - tele_location(1))>0.1) then
        print *, 'source depth is not consistent with the boundary data'
        pause
    end if
    
    read(10, *) bd_nr,bd_nt
    
    allocate(bd_time(bd_nr,bd_nt),bd_r(bd_nr),bd_t(bd_nt))
     
    read(10, *) bd_r
    read(10, *) bd_t

    do i=1,bd_nr
        do j=1,bd_nt
            read(10, *) bd_time(i,j)
        end do
    end do
    close(10)


    ! 统一计算震中距
    do i=1,nt
        do j=1,np
            lat1 = tt(i); lon1 = pp(j)
            lat2 = tele_location(2); lon2 = tele_location(3); 
            epicenter_dis(i,j)=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))/pi*180
        end do
    end do

    ! 指定深度 计算 走时
    idx0 = floor(((Radius-dep) -rr(1))/dr)+1
    r1 = min(1.0,((Radius-dep)-rr(idx0))/dr )
    L1_err=0; Linf_err=0

    do iit=1,nt
        do iip=1,np
            ! 震源位置: tele_location(1), 台站位置 Radius-dep，震中距：epicenter_dis(iit,iip),
            rec = Radius-dep
            dis = epicenter_dis(iit,iip)

            call interp2d(bd_r,bd_t,bd_nr,bd_nt,bd_time,rec,dis,traveltime)   

            error(iit,iip,1) = traveltime
            error(iit,iip,2) = (1-r1)*T(idx0,iit,iip)+r1*T(idx0+1,iit,iip)       
            error(iit,iip,3) = error(iit,iip,1) - error(iit,iip,2)
            L1_err = L1_err + abs(error(iit,iip,3))
            Linf_err = max(Linf_err, abs(error(iit,iip,3)))
        end do
    end do
    L1_err = L1_err/(nt*np)

    ! dep = 6071.0
    ! do iit=1,10
    !     call interp2d(bd_r,bd_t,bd_nr,bd_nt,bd_time,dep+iit,epicenter_dis(3,3),traveltime)  
    !     print *,  dep+iit,traveltime
    ! end do

    ! 释放动态数组
    deallocate(bd_time,bd_r,bd_t)
    
    write(*,'(a,f5.1,a,es10.3,a,es10.3,a)') 'L1 and Linf errors at depth ', dep, 'km are ', L1_err, ' s and ', Linf_err,' s.'

    


end subroutine
! ----------------------- end for ega3 ----------------------

! ----------------------- for ega4 ----------------------
subroutine eikonal_boundary_traveltime(rr,tt,pp,nr,nt,np,r0,t0,p0,velo_1d_fn,TTfield_path,isbd,bd_N,bd_S,bd_W,bd_E,bd_B,bd_T)
    ! input parameter
    integer :: nr,nt,np
    double precision :: rr(nr),tt(nt),pp(np)
    double precision :: r0,t0,p0
    character(len=80) :: velo_1d_fn,TTfield_path

    ! output parameter
    double precision :: bd_N(nr,np),bd_S(nr,np),bd_W(nr,nt),bd_E(nr,nt),bd_B(nt,np),bd_T(nt,np)
    logical :: isbd(5)

    ! 
    double precision :: degree(4),azimuthal(4)
    double precision :: maxdegree
    double precision :: epicenter_dis(nt,np) ! epicenter distance from tele source to study region

    ! mesh for 2d eikonal
    integer :: nr_2d,nt_2d
    double precision :: rr1_2d,rr2_2d,tt1_2d,tt2_2d,dr_2d,dt_2d,r0_2d,t0_2d
    double precision,allocatable :: rr_2d(:),tt_2d(:),spha_2d(:,:),sphb_2d(:,:)
    double precision,allocatable :: fun_2d(:,:),u_2d(:,:),T_2d(:,:)
    character(len=80) :: filename,form1,form2
    logical :: isfile,iscal
    integer :: file_nr_2d,file_nt_2d
    double precision,allocatable :: file_rr_2d(:),file_tt_2d(:)


    ! others
    double precision,parameter :: pi=3.14159265358979323846264338327950288
    double precision,parameter :: Earth_Radius = 6371.0
    double precision,parameter :: tol=0.001
    integer :: iir,iit,iip,iter  ! iteration index


    ! ---------- 首先，确定二维程函方程的计算范围 first, determine the region for 2d eikonal --------------------
    
    ! calculate the degree and azimuthal from the tele-source to the region
    ! S-W
    call Epicenter_Distance_sphere(tt(1),pp(1),t0,p0,degree(1))
    call Azimuthal_sphere(tt(1),pp(1),t0,p0,azimuthal(1))
    ! S-E
    call Epicenter_Distance_sphere(tt(1),pp(np),t0,p0,degree(2))
    call Azimuthal_sphere(tt(1),pp(np),t0,p0,azimuthal(2))
    ! N-W
    call Epicenter_Distance_sphere(tt(nt),pp(1),t0,p0,degree(3))
    call Azimuthal_sphere(tt(nt),pp(1),t0,p0,azimuthal(3))
    ! N-E
    call Epicenter_Distance_sphere(tt(nt),pp(np),t0,p0,degree(4))
    call Azimuthal_sphere(tt(nt),pp(np),t0,p0,azimuthal(4))
    
    ! print *, 'source (lat,lon) = ', t0*180/pi, ' ', p0*180/pi
    ! print *, 'degree is ', degree*180/pi
    ! print *, 'azimuthal is ', azimuthal*180/pi
    
    ! ---------- 然后，根据震中距和方位角，确定计算哪些边界，只有4种情况，出现其他情况，请报错 ------------
    ! ---------- determine to activate which boundaries. Only 4 cases, otherwise print error ------------
    maxdegree = max(degree(1),degree(2),degree(3),degree(4)) 
    if (maxdegree== degree(1)) then
        ! N + E
        isbd = [ .true. , .false. , .false. , .true. , .true. ] ! isbd: N,S,W,E,B
    else if (maxdegree == degree(2)) then
        ! N + W
        isbd = [ .true. , .false. , .true. , .false. , .true. ]
    else if (maxdegree == degree(3)) then
        ! S + E
        isbd = [ .false. , .true. , .false. , .true. , .true. ]
    else if (maxdegree == degree(4)) then
        ! S + W
        isbd = [ .false. , .true. , .true. , .false. , .true. ]
    end if

    ! ----------  if source is on the west or ease direction, N and S boundaries are deactivated
    if (t0 > tt(1) .and. t0 < tt(nt)) then
        isbd(1) = .false.
        isbd(2) = .false.
    end if
    ! ----------  if source is on the north or south direction, west and east boundaries are deactivated
    if (p0 > pp(1) .and. p0 < pp(np)) then
        isbd(3) = .false.
        isbd(4) = .false.
    end if


    ! ---------- 生成二维网格 ----------------------
    ! ---------- generate 2d mesh ----------------------
    dr_2d = 2.0; dt_2d = pi/1000; 
    rr1_2d = Earth_Radius-3000; rr2_2d = Earth_Radius+100;
    tt1_2d = -10.0/180.0*pi; tt2_2d = 100.0/180.0*pi;
    nr_2d = floor((rr2_2d-rr1_2d)/dr_2d)+1; 
    nt_2d = floor((tt2_2d-tt1_2d)/dt_2d)+1;

    allocate(rr_2d(nr_2d),tt_2d(nt_2d))

    do iir=1,nr_2d
        rr_2d(iir)=rr1_2d+(iir-1)*dr_2d
    end do
    do iit=1,nt_2d
        tt_2d(iit)=tt1_2d+(iit-1)*dt_2d
    end do

    

    ! ----------- 读取一维速度文件，产生速度场 -----------------------
    ! ----------- read the 1d velocity file -----------------------
    allocate(fun_2d(nr_2d,nt_2d),spha_2d(nr_2d,nt_2d),sphb_2d(nr_2d,nt_2d),u_2d(nr_2d,nt_2d),T_2d(nr_2d,nt_2d))
    
    call Read_1D_velocity_2d(velo_1d_fn,rr_2d,nr_2d,nt_2d,fun_2d)
    do iir=1,nr_2d
        do iit=1,nt_2d
            spha_2d(iir,iit)=1.0; 
            sphb_2d(iir,iit)=1.0;
            u_2d(iir,iit) = 0.0
            T_2d(iir,iit) = 0.0
        end do
    end do
    ! --------- 速度场文件输出 output 2d velocity file -------------
    filename=trim(TTfield_path)//'/model2d'
    inquire(file=filename, exist=isfile)
    if (.not. isfile) then
        open(10,file=filename)
        write(10,*) nr_2d,' ',nt_2d
        write(10,*) rr_2d
        write(10,*) tt_2d
        do iir =1,nr_2d
            do iit =1,nt_2d
                write(10,*) 1.0/fun_2d(iir,iit)
            end do
            print *, Earth_Radius - rr_2d(iir), 1.0/fun_2d(iir,1)
        end do
        close(10)
    end if
    
    ! ----------- 核心，计算对应远震产生的二维走时场 -------------------
    ! ----------- calculate the 2d traveltime field -----------------
    r0_2d = r0; t0_2d = 0.0
    write(*,'(a)') &
        & 'calculating 2d traveltime field, region: '
    write(*,'(a, f7.2,a,f7.2,a)') &
        & 'Depth: ', Earth_Radius-rr2_2d, ' km - ' ,Earth_Radius-rr1_2d, ' km' 
    write(*,'(a, f7.2,a,f7.2)') &
        & 'Epicenter dis: ', tt1_2d*180.0/pi, ' - ' ,tt2_2d*180.0/pi
    write(*,'(a,f7.2,f7.2,i5,i5)') &
        & '(dr,dt,nr,nt) = ', dr_2d,dt_2d*180.0/pi,nr_2d,nt_2d

    ! determine the 2d time file name
    write(form1,'(i4)') floor(r0)
    write(form2,'(i3)') int((r0-floor(r0))*1000)
    filename=trim(TTfield_path)//'/traveltime2d_'//trim(form1)//'_'//trim(form2)
    

    ! load or calculate the traveltime 
    inquire(file=filename, exist=isfile)
    if (isfile) then
        ! if the traveltime file exist, just read it
        write(*,'(a,a)') 'file exist, loading ... :', filename
        open(10,file=filename)
        read(10,*) file_nr_2d,file_nt_2d
        allocate(file_rr_2d(file_nr_2d),file_tt_2d(file_nt_2d))
        read(10,*) file_rr_2d
        read(10,*) file_tt_2d
        if (file_nr_2d == nr_2d .and. file_nt_2d == nt_2d .and. &
          & abs(file_rr_2d(1)-rr_2d(1))<tol .and. abs(file_tt_2d(1) - tt_2d(1))<tol .and. &
          & abs(file_rr_2d(file_nr_2d)-rr_2d(nr_2d))<tol .and. abs(file_tt_2d(file_nt_2d)-tt_2d(nt_2d))<tol) then
            do iir =1,nr_2d
                do iit =1,nt_2d
                    read(10,*) T_2d(iir,iit)
                end do
            end do
            write(*,'(a)') 'travletime filed loaded'
            iscal = .false.
        else
            ! otherwise, delete this file
            call system('rm '//trim(filename))
            write(*,'(a)') 'travletime filed error, need recalculation'
            iscal = .true.
        end if
        close(10)
        deallocate(file_rr_2d,file_tt_2d)

    else
        write(*,'(a)') 'file does not exist, calculating ... '
        iscal = .true.
    end if


    if (iscal) then
        ! calculate the traveltime field
        call FSM_WENO3_PS_sphe_2d(rr_2d,tt_2d,nr_2d,nt_2d,spha_2d,sphb_2d,T_2d,fun_2d,r0_2d,t0_2d,u_2d)

        ! then, generate the traveltime file
        open(10,file=filename)
        write(10,*) nr_2d,' ',nt_2d
        write(10,*) rr_2d
        write(10,*) tt_2d
        do iir =1,nr_2d
            do iit =1,nt_2d
                write(10,*) T_2d(iir,iit)
            end do
        end do
        close(10)
    end if
    ! ----------- now we calculate the 2d traveltime  T_2d(iir,iit)


    ! ----------- 插值到指定边界中 ------------------------
    ! ----------- interpolate T_2d onto the boundary condition bd_N - bd_B

    ! ----------- 计算震中距 calculate epicenter distance ----------------------
    do iit=1,nt
        do iip=1,np
            call Epicenter_Distance_sphere(tt(iit),pp(iip),t0,p0,epicenter_dis(iit,iip))
        end do
    end do

    ! interpolation for bd_N and bd_S
    do iir=1,nr
        do iip=1,np
            ! 南边界 South boundary
            call interp2d(rr_2d,tt_2d,nr_2d,nt_2d,T_2d,rr(iir),epicenter_dis(1,iip),bd_S(iir,iip))    
            ! 北边界 North boundary
            call interp2d(rr_2d,tt_2d,nr_2d,nt_2d,T_2d,rr(iir),epicenter_dis(nt,iip),bd_N(iir,iip)) 
        end do
    end do
    
    ! interpolation for bd_W and bd_E
    do iir = 1,nr
        do iit = 1,nt
            ! 西边界 West boundary
            call interp2d(rr_2d,tt_2d,nr_2d,nt_2d,T_2d,rr(iir),epicenter_dis(iit,1),bd_W(iir,iit))    
            ! 东边界 East boundary
            call interp2d(rr_2d,tt_2d,nr_2d,nt_2d,T_2d,rr(iir),epicenter_dis(iit,np),bd_E(iir,iit))      
        end do
    end do

    ! interpolation for bd_B
    do iit = 1,nt
        do iip = 1,np
            ! 底边界 bottom boundary
            call interp2d(rr_2d,tt_2d,nr_2d,nt_2d,T_2d,rr(1),epicenter_dis(iit,iip),bd_B(iit,iip))     
            ! surface boundary
            call interp2d(rr_2d,tt_2d,nr_2d,nt_2d,T_2d,Earth_Radius,epicenter_dis(iit,iip),bd_T(iit,iip))     
        end do
    end do

    ! ----------- release space --------------
    deallocate(rr_2d,tt_2d,spha_2d,sphb_2d,fun_2d,u_2d,T_2d)

end subroutine

subroutine FSM_WENO3_tele_sphe_3d_ver2(rr,tt,pp,nr,nt,np,spha,sphb,sphc,sphf,fun,isbd,bd_N,bd_S,bd_W,bd_E,bd_B,T)
    ! a -d -e
    ! -d b -f
    ! -e -f c
    integer nr,nt,np
    double precision :: dr,dt,dp,rr(nr),tt(nt),pp(np),a(nr,nt,np),b(nr,nt,np),c(nr,nt,np),f(nr,nt,np)
    double precision :: spha(nr,nt,np),sphb(nr,nt,np),sphc(nr,nt,np),sphf(nr,nt,np)
    double precision :: fun(nr,nt,np),T(nr,nt,np),ischange(nr,nt,np),xi(nr,nt,np),eta(nr,nt,np)
    double precision :: bd_N(nr,np),bd_S(nr,np),bd_W(nr,nt),bd_E(nr,nt),bd_B(nt,np)
    logical :: isbd(5)  ! N,S,W,E,B
    double precision :: T_old(nr,nt,np)
    double precision :: sigr,sigt,sigp,coe,pr1,pr2,wr1,wr2,pt1,pt2,wt1,wt2,pp1,pp2,wp1,wp2,tpT
    double precision :: L1_dif,Linf_dif,old_L1_dif
    double precision,parameter :: tol = (10.0)**(-5),eps=10.0**(-12)
    integer,parameter :: MaxIter=1000
    double precision,parameter :: infT = 2000.0
    integer iter,rdirec,tdirec,pdirec,iir,iit,iip,stop_count
    double precision :: tmp  
    
    stop_count = 0;

    ! ------------------------ 构造网格 ------------------------ 
    dr=rr(2)-rr(1); dt=tt(2)-tt(1); dp=pp(2)-pp(1)

    ! ------------------------ 构造矩阵 build eikonal matrix -------------------------  
    ! ------------------------ 初始化走时场 initialize traveltime field --------------
    !    a -d -e
    !   -d  b -f
    !   -e -f  c
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                a(iir,iit,iip) = spha(iir,iit,iip)
                b(iir,iit,iip) = sphb(iir,iit,iip)/(rr(iir)**2)
                c(iir,iit,iip) = sphc(iir,iit,iip)/(rr(iir)**2*cos(tt(iit))**2)
                f(iir,iit,iip) = sphf(iir,iit,iip)/(rr(iir)**2*cos(tt(iit)))
                ischange(iir,iit,iip) = 1
                T(iir,iit,iip) = infT
            end do
        end do
    end do

    ! ----------------------- 赋予边界值 -----------------------------
    ! N boundary
    if (isbd(1)) then
        do iir=1,nr
            do iip=1,np
                T(iir,nt,iip) = bd_N(iir,iip)
                ischange(iir,nt,iip) = 0
            end do
        end do
    end if
    ! S boundary
    if (isbd(2)) then
        do iir=1,nr
            do iip=1,np
                T(iir,1,iip) = bd_S(iir,iip)
                ischange(iir,1,iip) = 0
            end do
        end do
    end if
    ! W boundary
    if (isbd(3)) then
        do iir=1,nr
            do iit=1,nt
                T(iir,iit,1) = bd_W(iir,iit)
                ischange(iir,iit,1) = 0
            end do
        end do
    end if
    ! E boundary
    if (isbd(4)) then
        do iir=1,nr
            do iit=1,nt
                T(iir,iit,np) = bd_E(iir,iit)
                ischange(iir,iit,np) = 0
            end do
        end do
    end if
    ! Bottom boundary
    if (isbd(5)) then
        do iit=1,nt
            do iip=1,np
                T(1,iit,iip) = bd_B(iit,iip)
                ischange(1,iit,iip) = 0
            end do
        end do
    end if

    L1_dif = 0; old_L1_dif = 0;

    ! 正式迭代，更新 iteration start, update T
    do iter =1,MaxIter
        T_old = T
        if (L1_dif > old_L1_dif) then
            stop_count = stop_count+1
        end if
        old_L1_dif = L1_dif
        L1_dif=0; Linf_dif=0;
        do rdirec = -1,1,2
            do tdirec = -1,1,2
                do pdirec = -1,1,2

                    !x: nr <-> 1, y: nt <-> 1, z: np <-> 1

                    do iir=nint(0.5+nr/2.0+(nr/2.0-0.5)*rdirec),nint(0.5+nr/2.0+(-nr/2.0+0.5)*rdirec),-rdirec
                        do iit=nint(0.5+nt/2.0+(nt/2.0-0.5)*tdirec),nint(0.5+nt/2.0+(-nt/2.0+0.5)*tdirec),-tdirec
                            do iip=nint(0.5+np/2.0+(np/2.0-0.5)*pdirec),nint(0.5+np/2.0+(-np/2.0+0.5)*pdirec),-pdirec
                                

                                if(ischange(iir,iit,iip)==1) then
                                    if( iir>1 .and. iir<nr .and. iit>1 .and. iit<Nt .and. iip>1 .and. iip<Np) then

                                        sigr=1.05*sqrt(a(iir,iit,iip)); 
                                        sigt=1.05*sqrt(b(iir,iit,iip)); 
                                        sigp=1.05*sqrt(c(iir,iit,iip))
                                        coe=1.0/((sigr/dr)+(sigt/dt)+(sigp/dp))

                                        ! forward and backward partial derivatives

                                        if (iir==2) then
                                            pr1=(T(iir,iit,iip)-T(iir-1,iit,iip))/dr;
                                            wr2=1.0/(1.0+2*((eps+(T(iir,iit,iip)-2*T(iir+1,iit,iip)+T(iir+2,iit,iip))**2)/ &
                                                        & (eps+(T(iir-1,iit,iip)-2*T(iir,iit,iip)+T(iir+1,iit,iip))**2))**2)
                                            pr2=(1-wr2)*(T(iir+1,iit,iip)-T(iir-1,iit,iip))/2/dr+ &
                                                  & wr2*(-3*T(iir,iit,iip)+4*T(iir+1,iit,iip)-T(iir+2,iit,iip))/2/dr;
                                        elseif (iir==nr-1) then
                                            wr1=1.0/(1.0+2*((eps+(T(iir,iit,iip)-2*T(iir-1,iit,iip)+T(iir-2,iit,iip))**2)/ &
                                                        & (eps+(T(iir+1,iit,iip)-2*T(iir,iit,iip)+T(iir-1,iit,iip))**2))**2)
                                            pr1=(1-wr1)*(T(iir+1,iit,iip)-T(iir-1,iit,iip))/2/dr+ &
                                                  & wr1*(3*T(iir,iit,iip)-4*T(iir-1,iit,iip)+T(iir-2,iit,iip))/2/dr;
                                            pr2=(T(iir+1,iit,iip)-T(iir,iit,iip))/dr;                                       
                                        else
                                            wr1=1.0/(1.0+2*((eps+(T(iir,iit,iip)-2*T(iir-1,iit,iip)+T(iir-2,iit,iip))**2)/ &
                                                        & (eps+(T(iir+1,iit,iip)-2*T(iir,iit,iip)+T(iir-1,iit,iip))**2))**2)
                                            pr1=(1.0-wr1)*(T(iir+1,iit,iip)-T(iir-1,iit,iip))/2/dr+ &
                                                  & wr1*(3*T(iir,iit,iip)-4*T(iir-1,iit,iip)+T(iir-2,iit,iip))/2/dr;
                                            wr2=1.0/(1.0+2*((eps+(T(iir,iit,iip)-2*T(iir+1,iit,iip)+T(iir+2,iit,iip))**2)/ &
                                                        & (eps+(T(iir-1,iit,iip)-2*T(iir,iit,iip)+T(iir+1,iit,iip))**2))**2)
                                            pr2=(1.0-wr2)*(T(iir+1,iit,iip)-T(iir-1,iit,iip))/2/dr+ &
                                                  & wr2*(-3*T(iir,iit,iip)+4*T(iir+1,iit,iip)-T(iir+2,iit,iip))/2/dr;
                                        end if
    
                                        if (iit==2) then
                                            pt1=(T(iir,iit,iip)-T(iir,iit-1,iip))/dt; 
                                            wt2=1.0/(1.0+2*((eps+(T(iir,iit,iip)-2*T(iir,iit+1,iip)+T(iir,iit+2,iip))**2)/ &
                                                        & (eps+(T(iir,iit-1,iip)-2*T(iir,iit,iip)+T(iir,iit+1,iip))**2))**2)
                                            pt2=(1.0-wt2)*(T(iir,iit+1,iip)-T(iir,iit-1,iip))/2/dt+ &
                                                & wt2*(-3*T(iir,iit,iip)+4*T(iir,iit+1,iip)-T(iir,iit+2,iip))/2/dt; 
                                        elseif (iit==nt-1) then
                                            wt1=1.0/(1+2*((eps+(T(iir,iit,iip)-2*T(iir,iit-1,iip)+T(iir,iit-2,iip))**2)/ &
                                                        & (eps+(T(iir,iit+1,iip)-2*T(iir,iit,iip)+T(iir,iit-1,iip))**2))**2)
                                            pt1=(1.0-wt1)*(T(iir,iit+1,iip)-T(iir,iit-1,iip))/2/dt+ &
                                                & wt1*(3*T(iir,iit,iip)-4*T(iir,iit-1,iip)+T(iir,iit-2,iip))/2/dt;
                                            pt2=(T(iir,iit+1,iip)-T(iir,iit,iip))/dt;  
                                        else
                                            wt1=1.0/(1.0+2*((eps+(T(iir,iit,iip)-2*T(iir,iit-1,iip)+T(iir,iit-2,iip))**2)/ &
                                                        & (eps+(T(iir,iit+1,iip)-2*T(iir,iit,iip)+T(iir,iit-1,iip))**2))**2)
                                            pt1=(1.0-wt1)*(T(iir,iit+1,iip)-T(iir,iit-1,iip))/2/dt+ &
                                                & wt1*(3*T(iir,iit,iip)-4*T(iir,iit-1,iip)+T(iir,iit-2,iip))/2/dt;
                                            wt2=1.0/(1.0+2*((eps+(T(iir,iit,iip)-2*T(iir,iit+1,iip)+T(iir,iit+2,iip))**2)/ &
                                                        & (eps+(T(iir,iit-1,iip)-2*T(iir,iit,iip)+T(iir,iit+1,iip))**2))**2)
                                            pt2=(1-wt2)*(T(iir,iit+1,iip)-T(iir,iit-1,iip))/2/dt+ &
                                                & wt2*(-3*T(iir,iit,iip)+4*T(iir,iit+1,iip)-T(iir,iit+2,iip))/2/dt; 
                                        end if
    
                                        if (iip==2) then
                                            pp1=(T(iir,iit,iip)-T(iir,iit,iip-1))/dp;
                                            wp2=1.0/(1.0+2*((eps+(T(iir,iit,iip)-2*T(iir,iit,iip+1)+T(iir,iit,iip+2))**2)/ &
                                                        & (eps+(T(iir,iit,iip-1)-2*T(iir,iit,iip)+T(iir,iit,iip+1))**2))**2)
                                            pp2=(1.0-wp2)*(T(iir,iit,iip+1)-T(iir,iit,iip-1))/2/dp+ &
                                                & wp2*(-3*T(iir,iit,iip)+4*T(iir,iit,iip+1)-T(iir,iit,iip+2))/2/dp; 
                                        elseif (iip==np-1) then
                                            wp1=1.0/(1.0+2*((eps+(T(iir,iit,iip)-2*T(iir,iit,iip-1)+T(iir,iit,iip-2))**2)/ &
                                                        & (eps+(T(iir,iit,iip+1)-2*T(iir,iit,iip)+T(iir,iit,iip-1))**2))**2)
                                            pp1=(1.0-wp1)*(T(iir,iit,iip+1)-T(iir,iit,iip-1))/2/dp+ &
                                                & wp1*(3*T(iir,iit,iip)-4*T(iir,iit,iip-1)+T(iir,iit,iip-2))/2/dp;
                                            pp2=(T(iir,iit,iip+1)-T(iir,iit,iip))/dp; 
                                        else
                                            wp1=1.0/(1.0+2*((eps+(T(iir,iit,iip)-2*T(iir,iit,iip-1)+T(iir,iit,iip-2))**2)/ &
                                                        & (eps+(T(iir,iit,iip+1)-2*T(iir,iit,iip)+T(iir,iit,iip-1))**2))**2)
                                            pp1=(1.0-wp1)*(T(iir,iit,iip+1)-T(iir,iit,iip-1))/2/dp+ &
                                                & wp1*(3*T(iir,iit,iip)-4*T(iir,iit,iip-1)+T(iir,iit,iip-2))/2/dp;
                                            wp2=1.0/(1.0+2*((eps+(T(iir,iit,iip)-2*T(iir,iit,iip+1)+T(iir,iit,iip+2))**2)/ &
                                                        & (eps+(T(iir,iit,iip-1)-2*T(iir,iit,iip)+T(iir,iit,iip+1))**2))**2)
                                            pp2=(1.0-wp2)*(T(iir,iit,iip+1)-T(iir,iit,iip-1))/2/dp+ &
                                                & wp2*(-3*T(iir,iit,iip)+4*T(iir,iit,iip+1)-T(iir,iit,iip+2))/2/dp; 
                                        end if


                                        ! calculate  LF Hamiltonian 

                                        HT=sqrt( a(iir,iit,iip)*((pr1+pr2)/2)**2 + b(iir,iit,iip)*((pt1+pt2)/2)**2 &
                                        & + c(iir,iit,iip)*((pp1+pp2)/2)**2 -2*f(iir,iit,iip)*(pt1+pt2)/2*(pp1+pp2)/2 )
                                        
                                        ! 更新 update timetable
                                        tpT=coe*(fun(iir,iit,iip)-HT)  &
                                        & +coe*(sigr*(pr2-pr1)/2+sigt*(pt2-pt1)/2+sigp*(pp2-pp1)/2)+T(iir,iit,iip);

                                        !write(*,*) fun(iir,iit,iip),HT

                                    
                                        T(iir,iit,iip) = tpT
                                        

                                    elseif (iir==1) then
                                        T(1,iit,iip) = max(2*T(2,iit,iip)-T(3,iit,iip),T(3,iit,iip))
                                    elseif (iir==nr) then
                                        T(nr,iit,iip) = max(2*T(nr-1,iit,iip)-T(nr-2,iit,iip),T(nr-2,iit,iip))
                                    elseif (iit==1) then        
                                        T(iir,1,iip) = max(2*T(iir,2,iip)-T(iir,3,iip),T(iir,3,iip))
                                    elseif (iit==nt) then        
                                        T(iir,nt,iip) = max(2*T(iir,nt-1,iip)-T(iir,nt-2,iip),T(iir,nt-2,iip))
                                    elseif (iip==1) then        
                                        T(iir,iit,1) = max(2*T(iir,iit,2)-T(iir,iit,3),T(iir,iit,3))
                                    elseif (iip==Np) then 
                                        T(iir,iit,np) = max(2*T(iir,iit,np-1)-T(iir,iit,np-2),T(iir,iit,np-2))

                                    end if
                                end if

                                
                            end do
                        end do
                    end do    

                end do
            end do
        end do

        
        ! 统计误差，判断迭代终止条件
        do iir=1,nr
            do iit=1,nt
                do iip=1,np
                    L1_dif=L1_dif+abs(T(iir,iit,iip)-T_old(iir,iit,iip))
                    Linf_dif=max(Linf_dif,abs(T(iir,iit,iip)-T_old(iir,iit,iip)))
                end do               
            end do
        end do
        L1_dif=L1_dif/((nr)*(nt)*(np))


        ! ################  iteration information  #################
        if (abs(L1_dif)<tol) then
            write(*,*) 'iter ',iter,', T is steadt'
            exit
        elseif (stop_count >= 5) then
            write(*,*) 'iter ',iter,', T is oscillating'
            exit
        else
            ! write(*,*) 'iter ',iter,', T is changing, L1 dif = ', L1_dif,'L inf dif = ', Linf_dif
        end if

        if (iter==MaxIter) then    
            write(*,*) 'iter ',iter,', max iteration steps'
        end if

    end do

end subroutine

subroutine FSM_O1_Adj_tele_sphe_3d(rr,tt,pp,nr,nt,np,Table,TableADJ,zeta,xi,eta,Nrec,rrec,trec,prec,sourceADJ,isbd)

    integer :: nr,nt,np,Nrec
    double precision :: dr,dt,dp,rr(nr),tt(nt),pp(np)
    double precision :: Table(nr,nt,np),TableADJ(nr,nt,np),delta(nr,nt,np)
    double precision :: zeta(nr,nt,np),xi(nr,nt,np),eta(nr,nt,np)
    double precision :: rrec(Nrec),trec(Nrec),prec(Nrec),sourceADJ(Nrec)
    double precision :: a1,a1m(nr,nt,np),a1p(nr,nt,np),a2,a2m(nr,nt,np),a2p(nr,nt,np)
    double precision :: b1,b1m(nr,nt,np),b1p(nr,nt,np),b2,b2m(nr,nt,np),b2p(nr,nt,np)
    double precision :: c1,c1m(nr,nt,np),c1p(nr,nt,np),c2,c2m(nr,nt,np),c2p(nr,nt,np)
    double precision :: coe
    logical :: isbd(5)

    integer,parameter :: MaxIter=100
    double precision,parameter :: tol=10.0**(-6),eps=10.0**(-6)

    integer :: iir,iit,iip,idi,idj,idk,rdirec,tdirec,pdirec
    double precision :: r1,r2,r3,Linf_dif,tpTabldADJ,Hamilton
    double precision :: tmpr1,tmpr2,tmpt1,tmpt2

    ! #######  网格间距  mesh size  ########
    dr = rr(2)-rr(1) 
    dt = tt(2)-tt(1) 
    dp = pp(2)-pp(1)
    
    ! #################  构造源项 build the sourece term, the delta function  #################
    do iir = 1,nr
        do iit = 1,nt
            do iip= 1,np
                delta(iir,iit,iip) = 0
                TableADJ(iir,iit,iip) = 0
            end do
        end do
    end do




    ! ------- 遍历每个台站，给予delta函数贡献 loop each station to contribute the delta function -------
    do ir=1,Nrec  
        idi=floor((rrec(ir)-rr(1))/dr+1); 
        idj=floor((trec(ir)-tt(1))/dt+1);
        idk=floor((prec(ir)-pp(1))/dp+1);

        r1 = min(1.0,(rrec(ir)-rr(idi))/dr); 
        r2 = min(1.0,(trec(ir)-tt(idj))/dt); 
        r3 = min(1.0,(prec(ir)-pp(idk))/dp); 

        delta(idi,idj,idk) = delta(idi,idj,idk) + sourceADJ(ir)*(1-r1)*(1-r2)*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj,idk) = delta(idi+1,idj,idk) + sourceADJ(ir)*r1*(1-r2)*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi,idj+1,idk) = delta(idi,idj+1,idk) + sourceADJ(ir)*(1-r1)*r2*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj+1,idk) = delta(idi+1,idj+1,idk) + sourceADJ(ir)*r1*r2*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi,idj,idk+1) = delta(idi,idj,idk+1) + sourceADJ(ir)*(1-r1)*(1-r2)*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj,idk+1) = delta(idi+1,idj,idk+1) + sourceADJ(ir)*r1*(1-r2)*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi,idj+1,idk+1) = delta(idi,idj+1,idk+1) + sourceADJ(ir)*(1-r1)*r2*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj+1,idk+1) = delta(idi+1,idj+1,idk+1) + sourceADJ(ir)*r1*r2*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
    end do

    ! #################  构造方程系数 calculate coefficients of equations  #################
    do iir = 2,nr-1
        do iit = 2,nt-1
            do iip= 2,np-1
                tmpr1 = (rr(iir-1)+rr(iir))/2; tmpt1 = (tt(iit-1)+tt(iit))/2;
                tmpr2 = (rr(iir)+rr(iir+1))/2; tmpt2 = (tt(iit)+tt(iit-1))/2;

                a1 = -(1+zeta(iir-1,iit,iip)+zeta(iir,iit,iip))*(Table(iir,iit,iip)-Table(iir-1,iit,iip))/dr
                a1m(iir,iit,iip) = (a1-abs(a1))/2; a1p(iir,iit,iip) = (a1+abs(a1))/2;
                a2 = -(1+zeta(iir,iit,iip)+zeta(iir+1,iit,iip))*(Table(iir+1,iit,iip)-Table(iir,iit,iip))/dr
                a2m(iir,iit,iip) = (a2-abs(a2))/2; a2p(iir,iit,iip) = (a2+abs(a2))/2;

                !b1 = -(1-xi(iir,iit-1,iip)-xi(iir,iit,iip))/(tmpr1**2)*(Table(iir,iit,iip)-Table(iir,iit-1,iip))/dt &
                !   & -(eta(iir,iit-1,iip)+eta(iir,iit,iip))/(cos(tmpt1)*tmpr1**2)/(4*dp) &
                b1 = -(1-xi(iir,iit-1,iip)-xi(iir,iit,iip))/(rr(iir)**2)*(Table(iir,iit,iip)-Table(iir,iit-1,iip))/dt &
                   & -(eta(iir,iit-1,iip)+eta(iir,iit,iip))/(cos(tmpt1)*rr(iir)**2)/(4*dp) &
                   & *((Table(iir,iit-1,iip+1)-Table(iir,iit-1,iip-1))+(Table(iir,iit,iip+1)-Table(iir,iit,iip-1)))
                b1m(iir,iit,iip) = (b1-abs(b1))/2; b1p(iir,iit,iip) = (b1+abs(b1))/2;

                !b2 = -(1-xi(iir,iit,iip)-xi(iir,iit+1,iip))/(tmpr2**2)*(Table(iir,iit+1,iip)-Table(iir,iit,iip))/dt &
                !   & -(eta(iir,iit,iip)+eta(iir,iit+1,iip))/(cos(tmpt2)*tmpr2**2)/(4*dp) &
                b2 = -(1-xi(iir,iit,iip)-xi(iir,iit+1,iip))/(rr(iir)**2)*(Table(iir,iit+1,iip)-Table(iir,iit,iip))/dt &
                   & -(eta(iir,iit,iip)+eta(iir,iit+1,iip))/(cos(tmpt2)*rr(iir)**2)/(4*dp) &
                   & *((Table(iir,iit,iip+1)-Table(iir,iit,iip-1))+(Table(iir,iit+1,iip+1)-Table(iir,iit+1,iip-1)))
                b2m(iir,iit,iip) = (b2-abs(b2))/2; b2p(iir,iit,iip) = (b2+abs(b2))/2;

                !c1 = -(eta(iir,iit,iip-1)+eta(iir,iit,iip))/(cos(tmpt1)*tmpr1**2)/(4*dt) &
                c1 = -(eta(iir,iit,iip-1)+eta(iir,iit,iip))/(cos(tt(iit))*rr(iir)**2)/(4*dt) &
                   & *((Table(iir,iit+1,iip-1)-Table(iir,iit-1,iip-1))+(Table(iir,iit+1,iip)-Table(iir,iit-1,iip))) &
                !   & -(1+xi(iir,iit,iip-1)+xi(iir,iit,iip))/(cos(tmpt1)**2*tmpr1**2)*(Table(iir,iit,iip)-Table(iir,iit,iip-1))/dp
                   & -(1+xi(iir,iit,iip-1)+xi(iir,iit,iip))/ &
                   & (cos(tt(iit))**2*rr(iir)**2)*(Table(iir,iit,iip)-Table(iir,iit,iip-1))/dp
                c1m(iir,iit,iip) = (c1-abs(c1))/2; c1p(iir,iit,iip) = (c1+abs(c1))/2;

                !c2 = -(eta(iir,iit,iip)+eta(iir,iit,iip+1))/(cos(tmpt2)*tmpr2**2)/(4*dt) &
                c2 = -(eta(iir,iit,iip)+eta(iir,iit,iip+1))/(cos(tt(iit))*rr(iir)**2)/(4*dt) &
                   & *((Table(iir,iit+1,iip)-Table(iir,iit-1,iip))+(Table(iir,iit+1,iip+1)-Table(iir,iit-1,iip+1))) &
                !   & -(1+xi(iir,iit,iip)+xi(iir,iit,iip+1))/(cos(tmpt2)**2*tmpr2**2)*(Table(iir,iit,iip+1)-Table(iir,iit,iip))/dp
                   & -(1+xi(iir,iit,iip)+xi(iir,iit,iip+1)) &
                   & /(cos(tt(iit))**2*rr(iir)**2)*(Table(iir,iit,iip+1)-Table(iir,iit,iip))/dp
                c2m(iir,iit,iip) = (c2-abs(c2))/2; c2p(iir,iit,iip) = (c2+abs(c2))/2;
                
            end do
        end do
    end do

    ! ################ 使用FSM计算伴随场 calculate adjoint field by FSM ######################

    do iter =1,MaxIter
        Linf_dif=0;
        do rdirec = -1,1,2
            do tdirec = -1,1,2
                do pdirec = -1,1,2

                    !x: nr-1 <-> 2, y: nt-1 <-> 2, z: np-1 <-> 2
                    do iir=nint(0.5+nr/2.0+(nr/2.0-0.5)*rdirec),nint(0.5+nr/2.0+(-nr/2.0+0.5)*rdirec),-rdirec
                        do iit=nint(0.5+nt/2.0+(nt/2.0-0.5)*tdirec),nint(0.5+nt/2.0+(-nt/2.0+0.5)*tdirec),-tdirec
                            do iip=nint(0.5+np/2.0+(np/2.0-0.5)*pdirec),nint(0.5+np/2.0+(-np/2.0+0.5)*pdirec),-pdirec
                                
                                if (iir == 1) then
                                    ! bottom boundary
                                    if (isbd(5)) then
                                        if (TableADJ(3,iit,iip)>=0) then
                                            TableADJ(1,iit,iip) = max(0.0,2*TableADJ(2,iit,iip)-TableADJ(3,iit,iip))
                                        else
                                            TableADJ(1,iit,iip) = min(0.0,2*TableADJ(2,iit,iip)-TableADJ(3,iit,iip))
                                        end if
                                    else
                                        TableADJ(1,iit,iip) = 0.0
                                    end if

                                elseif (iir == nr) then
                                    ! Top boundary
                                    TableADJ(nr,iit,iip) = 0.0

                                elseif (iit == 1) then
                                    ! south boundary
                                    if (isbd(2)) then
                                        if (TableADJ(iir,3,iip)>=0) then
                                            TableADJ(iir,1,iip) = max(0.0,2*TableADJ(iir,2,iip)-TableADJ(iir,3,iip))
                                        else
                                            TableADJ(iir,1,iip) = min(0.0,2*TableADJ(iir,2,iip)-TableADJ(iir,3,iip))
                                        end if
                                    else
                                        TableADJ(iir,1,iip) = 0.0
                                    end if

                                elseif (iit == nt) then
                                    ! north boundary
                                    if (isbd(1)) then
                                        if (TableADJ(iir,nt-2,iip)>=0) then
                                            TableADJ(iir,nt,iip) = max(0.0,2*TableADJ(iir,nt-1,iip)-TableADJ(iir,nt-2,iip))
                                        else
                                            TableADJ(iir,nt,iip) = min(0.0,2*TableADJ(iir,nt-1,iip)-TableADJ(iir,nt-2,iip))
                                        end if
                                    else
                                        TableADJ(iir,nt,iip) = 0.0
                                    end if

                                elseif (iip == 1) then
                                    ! west boundary    
                                    if (isbd(3)) then
                                        if (TableADJ(iir,iit,3)>=0) then
                                            TableADJ(iir,iit,1) = max(0.0,2*TableADJ(iir,iit,2)-TableADJ(iir,iit,3))
                                        else
                                            TableADJ(iir,iit,1) = min(0.0,2*TableADJ(iir,iit,2)-TableADJ(iir,iit,3))
                                        end if
                                    else
                                        TableADJ(iir,iit,1) = 0.0
                                    end if
                                elseif (iip == np) then
                                    ! east boundary    
                                    if (isbd(4)) then
                                        if (TableADJ(iir,iit,np-2)>=0) then
                                            TableADJ(iir,iit,np) = max(0.0,2*TableADJ(iir,iit,np-1)-TableADJ(iir,iit,np-2))
                                        else
                                            TableADJ(iir,iit,np) = min(0.0,2*TableADJ(iir,iit,np-1)-TableADJ(iir,iit,np-2))
                                        end if
                                    else
                                        TableADJ(iir,iit,np) = 0.0
                                    end if
                                else
                                    coe = (a2p(iir,iit,iip)-a1m(iir,iit,iip))/dr &
                                    & +(b2p(iir,iit,iip)-b1m(iir,iit,iip))/dt &
                                    & +(c2p(iir,iit,iip)-c1m(iir,iit,iip))/dp
                                    
                                    if (abs(coe)<eps) then
                                        TableADJ(iir,iit,iip) = 0.0
                                    else
                                        Hamilton = &
                                        &   (TableADJ(iir-1,iit,iip)*a1p(iir,iit,iip)-TableADJ(iir+1,iit,iip)*a2m(iir,iit,iip))/dr &
                                        & + (TableADJ(iir,iit-1,iip)*b1p(iir,iit,iip)-TableADJ(iir,iit+1,iip)*b2m(iir,iit,iip))/dt &
                                        & + (TableADJ(iir,iit,iip-1)*c1p(iir,iit,iip)-TableADJ(iir,iit,iip+1)*c2m(iir,iit,iip))/dp

                                        tpTableADJ = (delta(iir,iit,iip)+Hamilton)/coe
                                        

                                        !if (abs(tpTableADJ-TableADJ(iir,iit,iip))>Linf_dif) then
                                        !    print *, iir,iit,iip,tpTableADJ,TableADJ(iir,iit,iip)
                                        !end if
                                        Linf_dif = max(Linf_dif, abs(tpTableADJ-TableADJ(iir,iit,iip)))
                                        
                                        TableADJ(iir,iit,iip) = tpTableADJ
                                    end if  
                                
                                end if

                            end do
                        end do
                    end do 

                end do
            end do
        end do

        

        if (abs(Linf_dif)<tol) then
            write(*,*) 'iter ',iter,', TableADJ is steady'
            exit
        else
            write(*,*) 'iter ',iter,', TableADJ is changing, continue ... '
        end if

        if (iter==MaxIter) then    
        !    write(*,*) 'iter ',iter,', max iteration steps'
        end if

    end do
    
end subroutine

subroutine Adjoint_Source_Ddt(Ttime,TtimeT,tst,pst,nst,max_deg,sourceADJ,obj)
    integer :: nst
    double precision :: Ttime(nst),TtimeT(nst),sourceADJ(nst),obj
    double precision :: tst(nst),pst(nst),max_deg,deg
    integer :: i,j
    ! double precision,parameter :: pi=3.141592653589793

    sourceADJ = 0
    obj = 0
    do i=1,nst-1
        do j=i+1,nst          
            call Epicenter_Distance_sphere(tst(i),pst(i),tst(j),pst(j),deg)                  
            if (deg < max_deg) then
                obj = obj + ((Ttime(i)-Ttime(j))-(TtimeT(i)-TtimeT(j)))**2
                sourceADJ(i) = sourceADJ(i)+((Ttime(i)-Ttime(j))-(TtimeT(i)-TtimeT(j)))
                sourceADJ(j) = sourceADJ(j)-((Ttime(i)-Ttime(j))-(TtimeT(i)-TtimeT(j)))
            end if
        end do
    end do

end subroutine

subroutine Epicenter_Distance_sphere(lat1,lon1,lat2,lon2,deg)
    double precision :: lat1,lon1,lat2,lon2,deg
    double precision,parameter :: pi=3.1415926535897932384626433

    if((lat1-lat2)**2+(lon1-lon2)**2<0.0001) then
        deg=0.0
    else
        deg = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))
    end if
end subroutine

subroutine Azimuthal_sphere(lat1,lon1,lat2,lon2,azi)
    double precision :: lat1,lon1,lat2,lon2,azi
    double precision,parameter :: pi=3.1415926535897932384626433
    double precision :: x,y
    if((lat1-lat2)**2+(lon1-lon2)**2<0.0001) then
        azi=0.0
    else
        y = sin(lon2 - lon1) * cos(lat2)
        x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon2 - lon1)
        azi = atan2(y, x) 
        if (azi < 0) then
            azi = azi+2*pi
        end if
    end if
end subroutine

subroutine read_model_3d(read_model_name,nr,nt,np,fun,xi,eta,zeta,a,b,c,f)
    integer :: nr,nt,np
    double precision :: fun(nr,nt,np),xi(nr,nt,np),eta(nr,nt,np),zeta(nr,nt,np)
    double precision :: a(nr,nt,np),b(nr,nt,np),c(nr,nt,np),f(nr,nt,np)
    character(len=80) :: read_model_name
    double precision :: dp_tmp
    integer :: i,j,k
    double precision,parameter :: gamma = 0.0 

    open(100,file=trim(read_model_name))
    do i=1,nr
        do j=1,nt
            do k=1,np
                read(100,*) dp_tmp,dp_tmp,dp_tmp,fun(i,j,k),xi(i,j,k),eta(i,j,k)
                zeta(i,j,k)=gamma*sqrt(eta(i,j,k)**2+xi(i,j,k)**2)
                a(i,j,k)=1.0+2*zeta(i,j,k);
                b(i,j,k)=1.0-2*xi(i,j,k);
                c(i,j,k)=1.0+2*xi(i,j,k);
                f(i,j,k)=-2*eta(i,j,k);
            end do
        end do
    end do
    close(100)

end subroutine

! multiplicative factors are more robust （串行版本，可用于震源并行）
subroutine FSM_WENO3_PS_sphe_3d_mul(rr,tt,pp,nr,nt,np,spha,sphb,sphc,sphf,T,fun,r0,t0,p0,u)
    ! a -d -e
    ! -d b -f
    ! -e -f c
    integer nr,nt,np
    double precision :: dr,dt,dp,rr(nr),tt(nt),pp(np),a(nr,nt,np),b(nr,nt,np),c(nr,nt,np),f(nr,nt,np)
    double precision :: spha(nr,nt,np),sphb(nr,nt,np),sphc(nr,nt,np),sphf(nr,nt,np)
    double precision :: fun(nr,nt,np),T(nr,nt,np),ischange(nr,nt,np),u(nr,nt,np),xi(nr,nt,np),eta(nr,nt,np)
    double precision :: r1,r2,r3,a0,b0,c0,d0,f0,T0v(nr,nt,np),T0r(nr,nt,np),T0t(nr,nt,np),T0p(nr,nt,np)
    double precision :: tau(nr,nt,np),tau_old(nr,nt,np),r0,t0,p0
    integer :: idr0,idt0,idp0
    double precision :: sigr,sigt,sigp,coe,pr1,pr2,wr1,wr2,pt1,pt2,wt1,wt2,pp1,pp2,wp1,wp2,tpT
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err
    double precision :: t1,p1
    double precision,parameter :: tol = (10.0)**(-5),eps=10.0**(-12)
    integer,parameter :: MaxIter=1000
    integer iter,rdirec,tdirec,pdirec,iir,iit,iip
    double precision :: x0,y0,z0,x,y,z,xst,yst,zst,e11,e12,e13,e21,e22,e23,e31,e32,e33
    double precision :: xstr,xstt,xstp,ystr,ystt,ystp,zstr,zstt,zstp,tmp

    ! ------------------------ 构造网格 ------------------------ 
    dr=rr(2)-rr(1); dt=tt(2)-tt(1); dp=pp(2)-pp(1)

    ! ------------------------ 构造矩阵 -------------------------  
    !    a -d -e
    !   -d  b -f
    !   -e -f  c

    ! ------------------------ 构造 T0 ------------------------ 

    ! 震源处参数离散化
    idr0=floor((r0-rr(1))/dr+1); idt0=floor((t0-tt(1))/dt+1); idp0=floor((p0-pp(1))/dp+1); 
    r1 = min(1.0,(r0-rr(idr0))/dr); r2 = min(1.0,(t0-tt(idt0))/dt); r3 = min(1.0,(p0-pp(idp0))/dp); 

    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                a(iir,iit,iip) = spha(iir,iit,iip)
                b(iir,iit,iip) = sphb(iir,iit,iip)/(rr(iir)**2)
                c(iir,iit,iip) = sphc(iir,iit,iip)/(rr(iir)**2*cos(tt(iit))**2)
                f(iir,iit,iip) = sphf(iir,iit,iip)/(rr(iir)**2*cos(tt(iit)))
            end do
        end do
    end do

    a0=(1-r1)*(1-r2)*(1-r3)*a(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*a(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*a(idr0,idt0+1,idp0)+(1-r1)*r2*r3*a(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*a(idr0+1,idt0,idp0)+r1*(1-r2)*r3*a(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*a(idr0+1,idt0+1,idp0)+r1*r2*r3*a(idr0+1,idt0+1,idp0+1)

    b0=(1-r1)*(1-r2)*(1-r3)*b(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*b(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*b(idr0,idt0+1,idp0)+(1-r1)*r2*r3*b(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*b(idr0+1,idt0,idp0)+r1*(1-r2)*r3*b(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*b(idr0+1,idt0+1,idp0)+r1*r2*r3*b(idr0+1,idt0+1,idp0+1)

    c0=(1-r1)*(1-r2)*(1-r3)*c(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*c(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*c(idr0,idt0+1,idp0)+(1-r1)*r2*r3*c(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*c(idr0+1,idt0,idp0)+r1*(1-r2)*r3*c(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*c(idr0+1,idt0+1,idp0)+r1*r2*r3*c(idr0+1,idt0+1,idp0+1)

    f0=(1-r1)*(1-r2)*(1-r3)*f(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*f(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*f(idr0,idt0+1,idp0)+(1-r1)*r2*r3*f(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*f(idr0+1,idt0,idp0)+r1*(1-r2)*r3*f(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*f(idr0+1,idt0+1,idp0)+r1*r2*r3*f(idr0+1,idt0+1,idp0+1)
    !d0=-d0

    fun0=(1-r1)*(1-r2)*(1-r3)*fun(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*fun(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*fun(idr0,idt0+1,idp0)+(1-r1)*r2*r3*fun(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*fun(idr0+1,idt0,idp0)+r1*(1-r2)*r3*fun(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*fun(idr0+1,idt0+1,idp0)+r1*r2*r3*fun(idr0+1,idt0+1,idp0+1)
    
    ! 构造T0
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                r1 = rr(iir)
                t1 = tt(iit)
                p1 = pp(iip)
       
                T0v(iir,iit,iip) = fun0*sqrt( 1.0/a0*(r1-r0)**2 + c0/(c0*b0-f0**2)*(t1-t0)**2 &
                                          & + b0/(c0*b0-f0**2)*(p1-p0)**2 + 2*f0/(c0*b0-f0**2)*(t1-t0)*(p1-p0) ) 
                ischange(iir,iit,iip)=1
                if ( T0v(iir,iit,iip)==0 ) then
                    T0r(iir,iit,iip) = 0
                    T0t(iir,iit,iip) = 0
                    T0p(iir,iit,iip) = 0
                    
                else
                    T0r(iir,iit,iip) = fun0**2*(1.0/a0*(r1-r0))/T0v(iir,iit,iip)
                    T0t(iir,iit,iip) = fun0**2*(c0/(c0*b0-f0**2)*(t1-t0)+f0/(c0*b0-f0**2)*(p1-p0))/T0v(iir,iit,iip)
                    T0p(iir,iit,iip) = fun0**2*(b0/(c0*b0-f0**2)*(p1-p0)+f0/(c0*b0-f0**2)*(t1-t0))/T0v(iir,iit,iip)    
                end if

                if ( abs((rr(iir)-r0)/dr)<=2 .and. abs((tt(iit)-t0)/dt)<=2 .and. abs((pp(iip)-p0)/dp)<= 2) then
                    tau(iir,iit,iip) = 1  !震源周围几个点，直接认为是常速度结构，给出解析解，即和T0相等
                    !tau(iir,iit,iip) = u(iir,iit,iip) - T0v(iir,iit,iip)
                    ischange(iir,iit,iip)=0
                    if (iir==1 .or. iir==nr .or. iit==1 .or. iit==nt .or. iip==1 .or. iip==np) then
                        write(*,*) 'source on the boundary, mesh error'
                        print *, rr(iir),r0,tt(iit),t0,pp(iip),p0
                        pause
                    end if
                    !write(*,*) iir-idr0,iit-idt0,iip-idp0,u(iir,iit,iip),T0v(iir,iit,iip),u(iir,iit,iip)-T0v(iir,iit,iip)
                else
                    tau(iir,iit,iip) = 1
                    ischange(iir,iit,iip)=1
                end if
            end do
        end do
    end do

    L1_err=0; Linf_err=0
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                L1_err=L1_err+abs(u(iir,iit,iip)-T0v(iir,iit,iip))
                Linf_err=max(Linf_err,abs(T0v(iir,iit,iip)-u(iir,iit,iip)))   
            end do            
        end do
    end do
    L1_err = L1_err/(nr*np*nt)
    ! print *, L1_err

    ! 正式迭代，更新tau
    do iter =1,MaxIter
        tau_old = tau
        L1_dif=0; Linf_dif=0;L1_err=0;Linf_err=0;
        do rdirec = -1,1,2
            do tdirec = -1,1,2
                do pdirec = -1,1,2

                    !x: nr-1 <-> 2, y: nt-1 <-> 2, z: np-1 <-> 2
                    do iir=nint(0.5+nr/2.0+(nr/2.0-1.5)*rdirec),nint(0.5+nr/2.0+(-nr/2.0+1.5)*rdirec),-rdirec
                        do iit=nint(0.5+nt/2.0+(nt/2.0-1.5)*tdirec),nint(0.5+nt/2.0+(-nt/2.0+1.5)*tdirec),-tdirec
                            do iip=nint(0.5+np/2.0+(np/2.0-1.5)*pdirec),nint(0.5+np/2.0+(-np/2.0+1.5)*pdirec),-pdirec
                                
                                if(ischange(iir,iit,iip)==1) then
                                    sigr = 1.0*sqrt(a(iir,iit,iip))*T0v(iir,iit,iip)
                                    sigt = 1.0*sqrt(b(iir,iit,iip))*T0v(iir,iit,iip)
                                    sigp = 1.0*sqrt(c(iir,iit,iip))*T0v(iir,iit,iip)
                                    coe=1.0/((sigr/dr)+(sigt/dt)+(sigp/dp))
                                    ! 构造单方向梯度 3阶 WENO 格式

                                    if (iir==2) then
                                        pr1=(tau(iir,iit,iip)-tau(iir-1,iit,iip))/dr;
                                        wr2=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir+1,iit,iip)+tau(iir+2,iit,iip))**2)/ &
                                                    & (eps+(tau(iir-1,iit,iip)-2*tau(iir,iit,iip)+tau(iir+1,iit,iip))**2))**2)
                                        pr2=(1-wr2)*(tau(iir+1,iit,iip)-tau(iir-1,iit,iip))/2/dr+ &
                                              & wr2*(-3*tau(iir,iit,iip)+4*tau(iir+1,iit,iip)-tau(iir+2,iit,iip))/2/dr;
                                    elseif (iir==nr-1) then
                                        wr1=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir-1,iit,iip)+tau(iir-2,iit,iip))**2)/ &
                                                    & (eps+(tau(iir+1,iit,iip)-2*tau(iir,iit,iip)+tau(iir-1,iit,iip))**2))**2)
                                        pr1=(1-wr1)*(tau(iir+1,iit,iip)-tau(iir-1,iit,iip))/2/dr+ &
                                              & wr1*(3*tau(iir,iit,iip)-4*tau(iir-1,iit,iip)+tau(iir-2,iit,iip))/2/dr;
                                        pr2=(tau(iir+1,iit,iip)-tau(iir,iit,iip))/dr;                                       
                                    else
                                        wr1=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir-1,iit,iip)+tau(iir-2,iit,iip))**2)/ &
                                                    & (eps+(tau(iir+1,iit,iip)-2*tau(iir,iit,iip)+tau(iir-1,iit,iip))**2))**2)
                                        pr1=(1.0-wr1)*(tau(iir+1,iit,iip)-tau(iir-1,iit,iip))/2/dr+ &
                                              & wr1*(3*tau(iir,iit,iip)-4*tau(iir-1,iit,iip)+tau(iir-2,iit,iip))/2/dr;
                                        wr2=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir+1,iit,iip)+tau(iir+2,iit,iip))**2)/ &
                                                    & (eps+(tau(iir-1,iit,iip)-2*tau(iir,iit,iip)+tau(iir+1,iit,iip))**2))**2)
                                        pr2=(1.0-wr2)*(tau(iir+1,iit,iip)-tau(iir-1,iit,iip))/2/dr+ &
                                              & wr2*(-3*tau(iir,iit,iip)+4*tau(iir+1,iit,iip)-tau(iir+2,iit,iip))/2/dr;
                                    end if

                                    if (iit==2) then
                                        pt1=(tau(iir,iit,iip)-tau(iir,iit-1,iip))/dt; 
                                        wt2=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit+1,iip)+tau(iir,iit+2,iip))**2)/ &
                                                    & (eps+(tau(iir,iit-1,iip)-2*tau(iir,iit,iip)+tau(iir,iit+1,iip))**2))**2)
                                        pt2=(1.0-wt2)*(tau(iir,iit+1,iip)-tau(iir,iit-1,iip))/2/dt+ &
                                            & wt2*(-3*tau(iir,iit,iip)+4*tau(iir,iit+1,iip)-tau(iir,iit+2,iip))/2/dt; 
                                    elseif (iit==nt-1) then
                                        wt1=1.0/(1+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit-1,iip)+tau(iir,iit-2,iip))**2)/ &
                                                    & (eps+(tau(iir,iit+1,iip)-2*tau(iir,iit,iip)+tau(iir,iit-1,iip))**2))**2)
                                        pt1=(1.0-wt1)*(tau(iir,iit+1,iip)-tau(iir,iit-1,iip))/2/dt+ &
                                            & wt1*(3*tau(iir,iit,iip)-4*tau(iir,iit-1,iip)+tau(iir,iit-2,iip))/2/dt;
                                        pt2=(tau(iir,iit+1,iip)-tau(iir,iit,iip))/dt;  
                                    else
                                        wt1=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit-1,iip)+tau(iir,iit-2,iip))**2)/ &
                                                    & (eps+(tau(iir,iit+1,iip)-2*tau(iir,iit,iip)+tau(iir,iit-1,iip))**2))**2)
                                        pt1=(1.0-wt1)*(tau(iir,iit+1,iip)-tau(iir,iit-1,iip))/2/dt+ &
                                            & wt1*(3*tau(iir,iit,iip)-4*tau(iir,iit-1,iip)+tau(iir,iit-2,iip))/2/dt;
                                        wt2=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit+1,iip)+tau(iir,iit+2,iip))**2)/ &
                                                    & (eps+(tau(iir,iit-1,iip)-2*tau(iir,iit,iip)+tau(iir,iit+1,iip))**2))**2)
                                        pt2=(1-wt2)*(tau(iir,iit+1,iip)-tau(iir,iit-1,iip))/2/dt+ &
                                            & wt2*(-3*tau(iir,iit,iip)+4*tau(iir,iit+1,iip)-tau(iir,iit+2,iip))/2/dt; 
                                    end if

                                    if (iip==2) then
                                        pp1=(tau(iir,iit,iip)-tau(iir,iit,iip-1))/dp;
                                        wp2=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit,iip+1)+tau(iir,iit,iip+2))**2)/ &
                                                    & (eps+(tau(iir,iit,iip-1)-2*tau(iir,iit,iip)+tau(iir,iit,iip+1))**2))**2)
                                        pp2=(1.0-wp2)*(tau(iir,iit,iip+1)-tau(iir,iit,iip-1))/2/dp+ &
                                            & wp2*(-3*tau(iir,iit,iip)+4*tau(iir,iit,iip+1)-tau(iir,iit,iip+2))/2/dp; 
                                    elseif (iip==np-1) then
                                        wp1=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit,iip-1)+tau(iir,iit,iip-2))**2)/ &
                                                    & (eps+(tau(iir,iit,iip+1)-2*tau(iir,iit,iip)+tau(iir,iit,iip-1))**2))**2)
                                        pp1=(1.0-wp1)*(tau(iir,iit,iip+1)-tau(iir,iit,iip-1))/2/dp+ &
                                            & wp1*(3*tau(iir,iit,iip)-4*tau(iir,iit,iip-1)+tau(iir,iit,iip-2))/2/dp;
                                        pp2=(tau(iir,iit,iip+1)-tau(iir,iit,iip))/dp; 
                                    else
                                        wp1=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit,iip-1)+tau(iir,iit,iip-2))**2)/ &
                                                    & (eps+(tau(iir,iit,iip+1)-2*tau(iir,iit,iip)+tau(iir,iit,iip-1))**2))**2)
                                        pp1=(1.0-wp1)*(tau(iir,iit,iip+1)-tau(iir,iit,iip-1))/2/dp+ &
                                            & wp1*(3*tau(iir,iit,iip)-4*tau(iir,iit,iip-1)+tau(iir,iit,iip-2))/2/dp;
                                        wp2=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit,iip+1)+tau(iir,iit,iip+2))**2)/ &
                                                    & (eps+(tau(iir,iit,iip-1)-2*tau(iir,iit,iip)+tau(iir,iit,iip+1))**2))**2)
                                        pp2=(1.0-wp2)*(tau(iir,iit,iip+1)-tau(iir,iit,iip-1))/2/dp+ &
                                            & wp2*(-3*tau(iir,iit,iip)+4*tau(iir,iit,iip+1)-tau(iir,iit,iip+2))/2/dp; 
                                    end if

                                    !计算 LF Hamiltonian 
                                    
                                    Htau = sqrt( &
                                    &   a(iir,iit,iip)*(T0r(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pr1+pr2)/2)**2 &
                                    & + b(iir,iit,iip)*(T0t(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pt1+pt2)/2)**2 &
                                    & + c(iir,iit,iip)*(T0p(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pp1+pp2)/2)**2 &
                                    &-2*f(iir,iit,iip)*(T0t(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pt1+pt2)/2) &
                                    &                 *(T0p(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pp1+pp2)/2) )
                                    
                                    if (isnan(Htau)) then
                                        print *, iir,iit,iip
                                        print *, a(iir,iit,iip)*(T0r(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pr1+pr2)/2)**2
                                        print *, b(iir,iit,iip)*(T0t(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pt1+pt2)/2)**2
                                        print *, c(iir,iit,iip)*(T0p(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pp1+pp2)/2)**2
                                        pause
                                    end if

                                    ! 更新 timetable
                                    tpT=coe*(fun(iir,iit,iip)-Htau)  &
                                     & +coe*(sigr*(pr2-pr1)/2+sigt*(pt2-pt1)/2+sigp*(pp2-pp1)/2)+tau(iir,iit,iip);

                                    !write(*,*) fun(iir,iit,iip),Htau

                                    tau(iir,iit,iip) = tpT
                                end if

                                
                            end do
                        end do
                    end do    

                    ! 处理边界

                    do iit=1,nt
                        do iip=1,np
                                tau(1,iit,iip) = max(2*tau(2,iit,iip)-tau(3,iit,iip),tau(3,iit,iip))
                                tau(nr,iit,iip) = max(2*tau(nr-1,iit,iip)-tau(nr-2,iit,iip),tau(nr-2,iit,iip))
                        end do
                    end do
                    do iir=1,nr
                        do iip=1,np
                                tau(iir,1,iip) = max(2*tau(iir,2,iip)-tau(iir,3,iip),tau(iir,3,iip))
                                tau(iir,nt,iip) = max(2*tau(iir,nt-1,iip)-tau(iir,nt-2,iip),tau(iir,nt-2,iip))
                        end do
                    end do
                    do iir=1,nr
                        do iit=1,nt
                                tau(iir,iit,1) = max(2*tau(iir,iit,2)-tau(iir,iit,3),tau(iir,iit,3))
                                tau(iir,iit,np) = max(2*tau(iir,iit,np-1)-tau(iir,iit,np-2),tau(iir,iit,np-2))
                        end do
                    end do   


                end do
            end do
        end do

        

        ! 统计误差，判断迭代终止条件
        do iir=2,nr-1
            do iit=2,nt-1
                do iip=2,np-1
                    L1_dif=L1_dif+abs(tau(iir,iit,iip)-tau_old(iir,iit,iip))*T0v(iir,iit,iip)
                    Linf_dif=max(Linf_dif,abs(tau(iir,iit,iip)-tau_old(iir,iit,iip))*T0v(iir,iit,iip))
                end do               
            end do
        end do
        
        do iir=2,nr-1
            do iit=2,nt-1
                do iip=2,np-1
                    L1_err=L1_err+abs(u(iir,iit,iip)-T0v(iir,iit,iip)*tau(iir,iit,iip))
                    Linf_err=max(Linf_err,abs(T0v(iir,iit,iip)*tau(iir,iit,iip)-u(iir,iit,iip)))   
                end do            
            end do
        end do
        L1_err=L1_err/((nr-2)*(nt-2)*(np-2))
        L1_dif=L1_dif/((nr-2)*(nt-2)*(np-2))

        if (abs(L1_dif)<tol ) then
            write(*,*) 'iter ',iter,', T is steadt'
            exit
        else
            ! write(*,*) 'iter ',iter,', T is changing, L1 dif = ', L1_dif,'L inf dif = ', Linf_dif
        end if

        if (iter==MaxIter) then    
            write(*,*) 'iter ',iter,', max iteration steps'
        end if
    end do

    T=T0v*tau
    
end subroutine

! 多重网格，更新参数。适用于反演网格深度非均匀的情况
subroutine Parameter_Update_Multigrid_ver2 &
        & (rr,tt,pp,nr,nt,np,all_Kernel,nk,invr,invt,invp,ninvr,ninvt,ninvp,ngrid,stepsize,update_value)
    ! forward modeling grid (fine grid)
    integer :: nr,nt,np,nk  ! number of kernels for all parameters
    double precision :: rr(nr),tt(nt),pp(np),dr,dt,dp
    double precision :: update_value(nr,nt,np,nk),para(nr,nt,np,nk),all_Kernel(nr,nt,np,nk)
    
    ! inversion grid (coarse grid)
    integer :: ninvr,ninvt,ninvp,ngrid   ! we have "ngrid" sets of inversion grid 
    double precision :: invr(ngrid,ninvr),invt(ngrid,ninvt),invp(ngrid,ninvp),dinvr,dinvt,dinvp
    double precision :: invkernel(ninvr,ninvt,ninvp,nk)
    double precision :: stepsize
    

    integer :: iir,iit,iip,iik
    integer :: idr,irt,idp,igrid
    double precision :: r1,r2,r3,Linf,pert

    dr=rr(2)-rr(1); dt=tt(2)-tt(1); dp=pp(2)-pp(1)

    
    

    ! 初始化 initiate kernel
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                do iik=1,nk
                    para(iir,iit,iip,iik) = 0
                    update_value(iir,iit,iip,iik) = 0
                end do
            end do
        end do
    end do
    
    

    
    do igrid=1,ngrid
        dinvt=invt(igrid,2)-invt(igrid,1); dinvp=invp(igrid,2)-invp(igrid,1)

        ! 初始化 initiate kernel
        do iir=1,ninvr
            do iit=1,ninvt
                do iip=1,ninvp
                    do iik=1,nk
                        invkernel(iir,iit,iip,iik) = 0
                    end do
                end do
            end do
        end do
        
        ! smooth the kernel
        do iir=1,nr

            r1 = -1
            do ii_invr = 1,ninvr-1
                if (rr(iir)> invr(igrid,ii_invr) .and. rr(iir)<=invr(igrid,ii_invr+1)) then
                    idr = ii_invr
                    r1 = (rr(iir) - invr(igrid,idr)) / (invr(igrid,idr+1) - invr(igrid,idr))  ! r1=0 -> r=invr(idr); r1=1 -> r=invr(idr+1)
                    exit
                end if
            end do
            if (r1 < 0) then
                cycle   ! 如果与反演网格点无关，考察下一个点 if r is out of inversion grid invr, turn to next point
            end if

            do iit=1,nt

                r2 = -1
                do ii_invt = 1,ninvt-1
                    if (tt(iit)> invt(igrid,ii_invt) .and. tt(iit)<=invt(igrid,ii_invt+1)) then
                        idt = ii_invt
                        r2 = (tt(iit) - invt(igrid,idt)) / (invt(igrid,idt+1) - invt(igrid,idt))  ! r2=0 -> t=invt(idt); r2=1 -> t=invt(idt+1)
                        exit
                    end if
                end do
                if (r2 < 0) then
                    cycle   ! 如果与反演网格点无关，考察下一个点 if t is out of inversion grid invt, turn to next point
                end if

                do iip=1,np
                    r3 = -1
                    do ii_invp = 1,ninvp-1
                        if (pp(iip)> invp(igrid,ii_invp) .and. pp(iip)<=invp(igrid,ii_invp+1)) then
                            idp = ii_invp
                            r3 = (pp(iip) - invp(igrid,idp)) / (invp(igrid,idp+1) - invp(igrid,idp))  ! r3=0 -> p=invp(idp); r3=1 -> p=invp(idp+1)
                            exit
                        end if
                    end do
                    if (r3 < 0) then
                        cycle   ! 如果与反演网格点无关，考察下一个点 if p is out of inversion grid invp, turn to next point
                    end if
                    
                    do iik=1,nk
                        invkernel(idr,idt,idp,iik) = invkernel(idr,idt,idp,iik) &
                                                & + (1-r1)*(1-r2)*(1-r3)*all_Kernel(iir,iit,iip,iik)
                        invkernel(idr+1,idt,idp,iik) = invkernel(idr+1,idt,idp,iik) &
                                                & + r1*(1-r2)*(1-r3)*all_Kernel(iir,iit,iip,iik)
                        invkernel(idr,idt+1,idp,iik) = invkernel(idr,idt+1,idp,iik) &
                                                & + (1-r1)*r2*(1-r3)*all_Kernel(iir,iit,iip,iik)
                        invkernel(idr+1,idt+1,idp,iik) = invkernel(idr+1,idt+1,idp,iik) &
                                                & + r1*r2*(1-r3)*all_Kernel(iir,iit,iip,iik)
                        invkernel(idr,idt,idp+1,iik) = invkernel(idr,idt,idp+1,iik) &
                                                & + (1-r1)*(1-r2)*r3*all_Kernel(iir,iit,iip,iik)
                        invkernel(idr+1,idt,idp+1,iik) = invkernel(idr+1,idt,idp+1,iik) &
                                                & + r1*(1-r2)*r3*all_Kernel(iir,iit,iip,iik)
                        invkernel(idr,idt+1,idp+1,iik) = invkernel(idr,idt+1,idp+1,iik) &
                                                & + (1-r1)*r2*r3*all_Kernel(iir,iit,iip,iik)
                        invkernel(idr+1,idt+1,idp+1,iik) = invkernel(idr+1,idt+1,idp+1,iik) &
                                                & + r1*r2*r3*all_Kernel(iir,iit,iip,iik)      
                    end do

                end do
            end do
        end do
        
        ! build update value
        do iir=1,nr
            r1 = -1
            do ii_invr = 1,ninvr-1
                if (rr(iir)> invr(igrid,ii_invr) .and. rr(iir)<=invr(igrid,ii_invr+1)) then
                    idr = ii_invr
                    r1 = (rr(iir) - invr(igrid,idr)) / (invr(igrid,idr+1) - invr(igrid,idr))  ! r1=0 -> r=invr(idr); r1=1 -> r=invr(idr+1)
                    exit
                end if
            end do
            if (r1 < 0) then
                cycle   ! 如果与反演网格点无关，考察下一个点 if r is out of inversion grid invr, turn to next point
            end if

            do iit=1,nt
                r2 = -1
                do ii_invt = 1,ninvt-1
                    if (tt(iit)> invt(igrid,ii_invt) .and. tt(iit)<=invt(igrid,ii_invt+1)) then
                        idt = ii_invt
                        r2 = (tt(iit) - invt(igrid,idt)) / (invt(igrid,idt+1) - invt(igrid,idt))  ! r2=0 -> t=invt(idt); r2=1 -> t=invt(idt+1)
                        exit
                    end if
                end do
                if (r2 < 0) then
                    cycle   ! 如果与反演网格点无关，考察下一个点 if t is out of inversion grid invt, turn to next point
                end if
        
                do iip=1,np
                    r3 = -1
                    do ii_invp = 1,ninvp-1
                        if (pp(iip)> invp(igrid,ii_invp) .and. pp(iip)<=invp(igrid,ii_invp+1)) then
                            idp = ii_invp
                            r3 = (pp(iip) - invp(igrid,idp)) / (invp(igrid,idp+1) - invp(igrid,idp))  ! r3=0 -> p=invp(idp); r3=1 -> p=invp(idp+1)
                            exit
                        end if
                    end do
                    if (r3 < 0) then
                        cycle   ! 如果与反演网格点无关，考察下一个点 if p is out of inversion grid invp, turn to next point
                    end if
                    
                    
                    do iik=1,nk
                        pert = 0.0
                        pert = pert + invkernel(idr,idt,idp,iik)*(1-r1)*(1-r2)*(1-r3) 
                        pert = pert + invkernel(idr+1,idt,idp,iik)*r1*(1-r2)*(1-r3) 
                        pert = pert + invkernel(idr,idt+1,idp,iik)*(1-r1)*r2*(1-r3) 
                        pert = pert + invkernel(idr+1,idt+1,idp,iik)*r1*r2*(1-r3) 
                        pert = pert + invkernel(idr,idt,idp+1,iik)*(1-r1)*(1-r2)*r3 
                        pert = pert + invkernel(idr+1,idt,idp+1,iik)*r1*(1-r2)*r3 
                        pert = pert + invkernel(idr,idt+1,idp+1,iik)*(1-r1)*r2*r3 
                        pert = pert + invkernel(idr+1,idt+1,idp+1,iik)*r1*r2*r3 
                        para(iir,iit,iip,iik) = para(iir,iit,iip,iik)+pert
                    end do
                    
                    
                end do
            end do
        end do


    end do

    
    ! rescale
    Linf = 0.0
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                do iik=1,nk
                    if (Linf < abs(para(iir,iit,iip,iik))) then
                        Linf = abs(para(iir,iit,iip,iik))
                    end if
                end do
            end do
        end do
    end do

    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                do iik=1,nk
                    update_value(iir,iit,iip,iik) = para(iir,iit,iip,iik)/Linf*stepsize
                end do
            end do
        end do
    end do

end subroutine

! ----------------------- end for ega4 ----------------------

! 并行版本
subroutine FSM_WENO3_PS_sphe_3d_mul_mpi(rr,tt,pp,nr,nt,np,spha,sphb,sphc,sphf,T,fun,r0,t0,p0,u)
    
    use mpi
    ! a -d -e
    ! -d b -f
    ! -e -f c
    integer nr,nt,np
    double precision :: dr,dt,dp,rr(nr),tt(nt),pp(np),a(nr,nt,np),b(nr,nt,np),c(nr,nt,np),f(nr,nt,np)
    double precision :: spha(nr,nt,np),sphb(nr,nt,np),sphc(nr,nt,np),sphf(nr,nt,np)
    double precision :: fun(nr,nt,np),T(nr,nt,np),ischange(nr,nt,np),u(nr,nt,np),xi(nr,nt,np),eta(nr,nt,np)
    double precision :: r1,r2,r3,a0,b0,c0,d0,f0,T0v(nr,nt,np),T0r(nr,nt,np),T0t(nr,nt,np),T0p(nr,nt,np)
    double precision :: tau(nr,nt,np),tau_old(nr,nt,np),r0,t0,p0
    integer :: idr0,idt0,idp0
    double precision :: sigr,sigt,sigp,coe,pr1,pr2,wr1,wr2,pt1,pt2,wt1,wt2,pp1,pp2,wp1,wp2,tpT
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err
    double precision :: t1,p1
    double precision,parameter :: tol = (10.0)**(-5),eps=10.0**(-12)
    integer,parameter :: MaxIter=1000
    integer iter,rdirec,tdirec,pdirec,iir,iit,iip
    double precision :: x0,y0,z0,x,y,z,xst,yst,zst,e11,e12,e13,e21,e22,e23,e31,e32,e33
    double precision :: xstr,xstt,xstp,ystr,ystt,ystp,zstr,zstt,zstp
    integer :: ileft,iright,jleft,jright,ii,jj,kk
    logical :: isexit

    integer :: ierr,myid,nproc,tag,istat(mpi_status_size),iproc,ptid(nr,np,nt)
    integer :: mpi_ptijk(3,np*nt),mpi_ptn,ipt,int_temp,my_npt
    double precision :: mpi_ptv(np*nt),dp_temp

    integer :: ptr1,ptr2,iNumber,total_ptn,Nptn,sNptn,tpNptn
    double precision :: ave_ni
    integer,allocatable :: proc_irange(:,:,:)


    isexit = .false.
    tag=99
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    
    ! ------------------------ 构造网格 ------------------------ 
    dr=rr(2)-rr(1); dt=tt(2)-tt(1); dp=pp(2)-pp(1)

    ! ------------------------ 构造矩阵 -------------------------  
    !    a -d -e
    !   -d  b -f
    !   -e -f  c

    ! ------------------------ 构造 T0 ------------------------ 
    
    ! 震源处参数离散化
    idr0=floor((r0-rr(1))/dr+1); idt0=floor((t0-tt(1))/dt+1); idp0=floor((p0-pp(1))/dp+1); 
    r1 = min(1.0,(r0-rr(idr0))/dr); r2 = min(1.0,(t0-tt(idt0))/dt); r3 = min(1.0,(p0-pp(idp0))/dp); 

    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                a(iir,iit,iip) = spha(iir,iit,iip)
                b(iir,iit,iip) = sphb(iir,iit,iip)/(rr(iir)**2)
                c(iir,iit,iip) = sphc(iir,iit,iip)/(rr(iir)**2*cos(tt(iit))**2)
                f(iir,iit,iip) = sphf(iir,iit,iip)/(rr(iir)**2*cos(tt(iit)))
            end do
        end do
    end do
    
    a0=(1-r1)*(1-r2)*(1-r3)*a(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*a(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*a(idr0,idt0+1,idp0)+(1-r1)*r2*r3*a(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*a(idr0+1,idt0,idp0)+r1*(1-r2)*r3*a(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*a(idr0+1,idt0+1,idp0)+r1*r2*r3*a(idr0+1,idt0+1,idp0+1)

    b0=(1-r1)*(1-r2)*(1-r3)*b(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*b(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*b(idr0,idt0+1,idp0)+(1-r1)*r2*r3*b(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*b(idr0+1,idt0,idp0)+r1*(1-r2)*r3*b(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*b(idr0+1,idt0+1,idp0)+r1*r2*r3*b(idr0+1,idt0+1,idp0+1)

    c0=(1-r1)*(1-r2)*(1-r3)*c(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*c(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*c(idr0,idt0+1,idp0)+(1-r1)*r2*r3*c(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*c(idr0+1,idt0,idp0)+r1*(1-r2)*r3*c(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*c(idr0+1,idt0+1,idp0)+r1*r2*r3*c(idr0+1,idt0+1,idp0+1)

    f0=(1-r1)*(1-r2)*(1-r3)*f(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*f(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*f(idr0,idt0+1,idp0)+(1-r1)*r2*r3*f(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*f(idr0+1,idt0,idp0)+r1*(1-r2)*r3*f(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*f(idr0+1,idt0+1,idp0)+r1*r2*r3*f(idr0+1,idt0+1,idp0+1)
    !d0=-d0

    fun0=(1-r1)*(1-r2)*(1-r3)*fun(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*fun(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*fun(idr0,idt0+1,idp0)+(1-r1)*r2*r3*fun(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*fun(idr0+1,idt0,idp0)+r1*(1-r2)*r3*fun(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*fun(idr0+1,idt0+1,idp0)+r1*r2*r3*fun(idr0+1,idt0+1,idp0+1)

    

    ! 构造T0
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                r1 = rr(iir)
                t1 = tt(iit)
                p1 = pp(iip)
       
                T0v(iir,iit,iip) = fun0*sqrt( 1.0/a0*(r1-r0)**2 + c0/(c0*b0-f0**2)*(t1-t0)**2 &
                                          & + b0/(c0*b0-f0**2)*(p1-p0)**2 + 2*f0/(c0*b0-f0**2)*(t1-t0)*(p1-p0) ) 
                ischange(iir,iit,iip)=1
                if ( T0v(iir,iit,iip)==0 ) then
                    T0r(iir,iit,iip) = 0
                    T0t(iir,iit,iip) = 0
                    T0p(iir,iit,iip) = 0
                else
                    T0r(iir,iit,iip) = fun0**2*(1.0/a0*(r1-r0))/T0v(iir,iit,iip)
                    T0t(iir,iit,iip) = fun0**2*(c0/(c0*b0-f0**2)*(t1-t0)+f0/(c0*b0-f0**2)*(p1-p0))/T0v(iir,iit,iip)
                    T0p(iir,iit,iip) = fun0**2*(b0/(c0*b0-f0**2)*(p1-p0)+f0/(c0*b0-f0**2)*(t1-t0))/T0v(iir,iit,iip)    
                end if


                if ( abs((rr(iir)-r0)/dr)<=2 .and. abs((tt(iit)-t0)/dt)<=2 .and. abs((pp(iip)-p0)/dp)<= 2) then
                    tau(iir,iit,iip) = 1  !震源周围几个点，直接认为是常速度结构，给出解析解，即和T0相等
                    !tau(iir,iit,iip) = u(iir,iit,iip) - T0v(iir,iit,iip)
                    ischange(iir,iit,iip)=0
                    if (iir==1 .or. iir==nr .or. iit==1 .or. iit==nt .or. iip==1 .or. iip==np) then
                        write(*,*) 'source on the boundary, mesh error'
                        print *, rr(iir),iir,t0*180/3.1415927,iit,p0*180/3.1415927,iip
                    !    pause
                    end if
                    !write(*,*) iir-idr0,iit-idt0,iip-idp0,u(iir,iit,iip),T0v(iir,iit,iip),u(iir,iit,iip)-T0v(iir,iit,iip)
                else
                    !tau(iir,iit,iip) = T(iir,iit,iip)/T0v(iir,iit,iip)
                    tau(iir,iit,iip) = 1
                    ischange(iir,iit,iip)=1
                end if
            end do
        end do
    end do

    L1_err=0; Linf_err=0
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                L1_err=L1_err+abs(u(iir,iit,iip)-T0v(iir,iit,iip))
                Linf_err=max(Linf_err,abs(T0v(iir,iit,iip)-u(iir,iit,iip)))   
                if(isnan(L1_err)) then
                    print *, u(iir,iit,iip),T0v(iir,iit,iip)
                    pause
                end if
            end do            
        end do
    end do
    L1_err = L1_err/(nr*np*nt)
    print *, L1_err


    allocate(proc_irange(2,nproc,nr+nt+np-3-6+1))

    ! 给每个线程分配计算点
    Nptn=0
    sNptn=0

    do level = 6,nr+nt+np-3
        ileft = max(2,level-np-nt+2)
        iright = min(level-4,nr-1)

        iNumber = iright - ileft + 1

        total_ptn = 0
        do ip = ileft,iright
            total_ptn = total_ptn+(min(level-ip-2,nt-1)-max(2,level-ip-np+1)+1)
        end do
        Nptn = Nptn + total_ptn  ! 总点数
        ave_ni = total_ptn*1.0/nproc

        if (iNumber <= nproc) then
            ! 即i的数量小于等于线程数，每个线程分配一个i就好
            tpNptn = 0
            MaxtpNptn = 0
            do ip = ileft,iright
                proc_irange(:,ip-ileft+1,level-5)=(/ip,ip/)
                tpNptn = (min(level-ip-2,nt-1)-max(2,level-ip-np+1)+1)
                MaxtpNptn = max(MaxtpNptn,tpNptn) 
            end do
            do ip = iright+1,ileft+nproc-1
                proc_irange(:,ip-ileft+1,level-5)=(/-1,-2/)
            end do
            sNptn = sNptn + MaxtpNptn
        else
            !if (myid .eq. 0) then
            !    print *, ileft,iright
            !end if
            !call sleep(2)
            
            ! i的数量大于线程数，平均分配，接近 total_ptn/nproc
            int_temp = 0
            tpNptn = 0
            MaxtpNptn = 0
            ptr1 = ileft
            iproc = 1
            do ptr2 = ileft,iright             
                int_temp = int_temp + (min(level-ptr2-2,nt-1)-max(2,level-ptr2-np+1)+1)
                tpNptn = tpNptn + (min(level-ptr2-2,nt-1)-max(2,level-ptr2-np+1)+1)
                !if (myid .eq. 0) then
                !    print *, ptr2,(min(level-ptr2-2,nt-1)-max(2,level-ptr2-np+1)+1),int_temp,ave_ni
                !end if
                !call sleep(2)
                if ((int_temp>=ave_ni*iproc) .or. (ptr2 .eq.iright)) then
                    !if (iproc>nproc) then
                    !    print *, iproc,nproc
                    !end if
                    proc_irange(:,iproc,level-5) = (/ptr1,ptr2/)
                    ptr1 = ptr2+1
                    iproc = iproc +1
                    MaxtpNptn = max(MaxtpNptn,tpNptn) 
                    tpNptn = 0
                end if
            end do   
            sNptn = sNptn + MaxtpNptn         
            !call sleep(2)
        end if
        !if (myid .eq. 0) then
        !    print *, proc_irange(1,:,level-5)
        !    print *, proc_irange(2,:,level-5)
        !    print *, ' '
        !end if
    end do

    ! print *, Nptn*1.0/sNptn


    ! 正式迭代，更新tau
    do iter =1,MaxIter
        tau_old = tau
        L1_dif=0; Linf_dif=0;L1_err=0;Linf_err=0;
        
        do rdirec = -1,1,2
            do tdirec = -1,1,2
                do pdirec = -1,1,2
                    my_npt = 0
                    do level = 6,nr+nt+np-3
                        ! 2<= ir <= nr-1;   2<= it <= nr-1;   2<=ip <=np-1

                        ileft = max(2,level-np-nt+2)
                        iright = min(level-4,nr-1)

                        mpi_ptn=0
                        do ii = proc_irange(1,myid+1,level-5),proc_irange(2,myid+1,level-5)
                        
                            jleft = max(2,level-ii-np+1)
                            jright = min(level-ii-2,nt-1)
                            do jj = jleft,jright

                                kk = level-ii-jj
                                
                                if(rdirec<0) then
                                    iir=ii
                                else
                                    iir=nr+1-ii
                                end if
                                if(tdirec<0) then
                                    iit=jj
                                else
                                    iit=nt+1-jj
                                end if
                                if(pdirec<0) then
                                    iip=kk
                                else
                                    iip=np+1-kk
                                end if
            

                                if(ischange(iir,iit,iip)==1) then
                                    

                                    !sigr=2*sqrt(a(iir,iit,iip)); sigt=2*sqrt(b(iir,iit,iip)); sigp=2*sqrt(c(iir,iit,iip))
                                    !coe=1.0/((sigr/dr)+(sigt/dt)+(sigp/dp))

                                    sigr = 1.0*sqrt(a(iir,iit,iip))*T0v(iir,iit,iip)
                                    sigt = 1.0*sqrt(b(iir,iit,iip))*T0v(iir,iit,iip)
                                    sigp = 1.0*sqrt(c(iir,iit,iip))*T0v(iir,iit,iip)
                                    coe=1.0/((sigr/dr)+(sigt/dt)+(sigp/dp))

                                    ! 构造单方向梯度 3阶 WENO 格式

                                    if (iir==2) then
                                        pr1=(tau(iir,iit,iip)-tau(iir-1,iit,iip))/dr;
                                        wr2=1.0/(1+2*((eps+(tau(iir,iit,iip)-2*tau(iir+1,iit,iip)+tau(iir+2,iit,iip))**2)/ &
                                                    & (eps+(tau(iir-1,iit,iip)-2*tau(iir,iit,iip)+tau(iir+1,iit,iip))**2))**2)
                                        pr2=(1-wr2)*(tau(iir+1,iit,iip)-tau(iir-1,iit,iip))/2/dr+ &
                                            & wr2*(-3*tau(iir,iit,iip)+4*tau(iir+1,iit,iip)-tau(iir+2,iit,iip))/2/dr;
                                    elseif (iir==nr-1) then
                                        wr1=1.0/(1+2*((eps+(tau(iir,iit,iip)-2*tau(iir-1,iit,iip)+tau(iir-2,iit,iip))**2)/ &
                                                    & (eps+(tau(iir+1,iit,iip)-2*tau(iir,iit,iip)+tau(iir-1,iit,iip))**2))**2)
                                        pr1=(1-wr1)*(tau(iir+1,iit,iip)-tau(iir-1,iit,iip))/2/dr+ &
                                            & wr1*(3*tau(iir,iit,iip)-4*tau(iir-1,iit,iip)+tau(iir-2,iit,iip))/2/dr;
                                        pr2=(tau(iir+1,iit,iip)-tau(iir,iit,iip))/dr;                                       
                                    else
                                        wr1=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir-1,iit,iip)+tau(iir-2,iit,iip))**2)/ &
                                                    & (eps+(tau(iir+1,iit,iip)-2*tau(iir,iit,iip)+tau(iir-1,iit,iip))**2))**2)
                                        pr1=(1.0-wr1)*(tau(iir+1,iit,iip)-tau(iir-1,iit,iip))/2/dr+ &
                                            & wr1*(3*tau(iir,iit,iip)-4*tau(iir-1,iit,iip)+tau(iir-2,iit,iip))/2/dr;
                                        wr2=1.0/(1.0+2*((eps+(tau(iir,iit,iip)-2*tau(iir+1,iit,iip)+tau(iir+2,iit,iip))**2)/ &
                                                    & (eps+(tau(iir-1,iit,iip)-2*tau(iir,iit,iip)+tau(iir+1,iit,iip))**2))**2)
                                        pr2=(1.0-wr2)*(tau(iir+1,iit,iip)-tau(iir-1,iit,iip))/2/dr+ &
                                            & wr2*(-3*tau(iir,iit,iip)+4*tau(iir+1,iit,iip)-tau(iir+2,iit,iip))/2/dr;
                                    end if

                                    if (iit==2) then
                                        pt1=(tau(iir,iit,iip)-tau(iir,iit-1,iip))/dt; 
                                        wt2=1.0/(1+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit+1,iip)+tau(iir,iit+2,iip))**2)/ &
                                                    & (eps+(tau(iir,iit-1,iip)-2*tau(iir,iit,iip)+tau(iir,iit+1,iip))**2))**2)
                                        pt2=(1-wt2)*(tau(iir,iit+1,iip)-tau(iir,iit-1,iip))/2/dt+ &
                                            & wt2*(-3*tau(iir,iit,iip)+4*tau(iir,iit+1,iip)-tau(iir,iit+2,iip))/2/dt; 
                                    elseif (iit==nt-1) then
                                        wt1=1.0/(1+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit-1,iip)+tau(iir,iit-2,iip))**2)/ &
                                                    & (eps+(tau(iir,iit+1,iip)-2*tau(iir,iit,iip)+tau(iir,iit-1,iip))**2))**2)
                                        pt1=(1-wt1)*(tau(iir,iit+1,iip)-tau(iir,iit-1,iip))/2/dt+ &
                                            & wt1*(3*tau(iir,iit,iip)-4*tau(iir,iit-1,iip)+tau(iir,iit-2,iip))/2/dt;
                                        pt2=(tau(iir,iit+1,iip)-tau(iir,iit,iip))/dt;  
                                    else
                                        wt1=1.0/(1+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit-1,iip)+tau(iir,iit-2,iip))**2)/ &
                                                    & (eps+(tau(iir,iit+1,iip)-2*tau(iir,iit,iip)+tau(iir,iit-1,iip))**2))**2)
                                        pt1=(1-wt1)*(tau(iir,iit+1,iip)-tau(iir,iit-1,iip))/2/dt+ &
                                            & wt1*(3*tau(iir,iit,iip)-4*tau(iir,iit-1,iip)+tau(iir,iit-2,iip))/2/dt;
                                        wt2=1.0/(1+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit+1,iip)+tau(iir,iit+2,iip))**2)/ &
                                                    & (eps+(tau(iir,iit-1,iip)-2*tau(iir,iit,iip)+tau(iir,iit+1,iip))**2))**2)
                                        pt2=(1-wt2)*(tau(iir,iit+1,iip)-tau(iir,iit-1,iip))/2/dt+ &
                                            & wt2*(-3*tau(iir,iit,iip)+4*tau(iir,iit+1,iip)-tau(iir,iit+2,iip))/2/dt; 
                                    end if

                                    if (iip==2) then
                                        pp1=(tau(iir,iit,iip)-tau(iir,iit,iip-1))/dp;
                                        wp2=1.0/(1+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit,iip+1)+tau(iir,iit,iip+2))**2)/ &
                                                    & (eps+(tau(iir,iit,iip-1)-2*tau(iir,iit,iip)+tau(iir,iit,iip+1))**2))**2)
                                        pp2=(1-wp2)*(tau(iir,iit,iip+1)-tau(iir,iit,iip-1))/2/dp+ &
                                            & wp2*(-3*tau(iir,iit,iip)+4*tau(iir,iit,iip+1)-tau(iir,iit,iip+2))/2/dp; 
                                    elseif (iip==np-1) then
                                        wp1=1.0/(1+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit,iip-1)+tau(iir,iit,iip-2))**2)/ &
                                                    & (eps+(tau(iir,iit,iip+1)-2*tau(iir,iit,iip)+tau(iir,iit,iip-1))**2))**2)
                                        pp1=(1-wp1)*(tau(iir,iit,iip+1)-tau(iir,iit,iip-1))/2/dp+ &
                                            & wp1*(3*tau(iir,iit,iip)-4*tau(iir,iit,iip-1)+tau(iir,iit,iip-2))/2/dp;
                                        pp2=(tau(iir,iit,iip+1)-tau(iir,iit,iip))/dp; 
                                    else
                                        wp1=1.0/(1+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit,iip-1)+tau(iir,iit,iip-2))**2)/ &
                                                    & (eps+(tau(iir,iit,iip+1)-2*tau(iir,iit,iip)+tau(iir,iit,iip-1))**2))**2)
                                        pp1=(1-wp1)*(tau(iir,iit,iip+1)-tau(iir,iit,iip-1))/2/dp+ &
                                            & wp1*(3*tau(iir,iit,iip)-4*tau(iir,iit,iip-1)+tau(iir,iit,iip-2))/2/dp;
                                        wp2=1.0/(1+2*((eps+(tau(iir,iit,iip)-2*tau(iir,iit,iip+1)+tau(iir,iit,iip+2))**2)/ &
                                                    & (eps+(tau(iir,iit,iip-1)-2*tau(iir,iit,iip)+tau(iir,iit,iip+1))**2))**2)
                                        pp2=(1-wp2)*(tau(iir,iit,iip+1)-tau(iir,iit,iip-1))/2/dp+ &
                                            & wp2*(-3*tau(iir,iit,iip)+4*tau(iir,iit,iip+1)-tau(iir,iit,iip+2))/2/dp; 
                                    end if

                                    !计算 LF Hamiltonian 

                                    !Htau=sqrt( a(iir,iit,iip)*((pr1+pr2)/2)**2 + b(iir,iit,iip)*((pt1+pt2)/2)**2 &
                                    !& + c(iir,iit,iip)*((pp1+pp2)/2)**2 -2*f(iir,iit,iip)*(pt1+pt2)/2*(pp1+pp2)/2 &
                                    !& + a(iir,iit,iip)*(pr1+pr2)*T0r(iir,iit,iip) &
                                    !& + b(iir,iit,iip)*(pt1+pt2)*T0t(iir,iit,iip) - f(iir,iit,iip)*(pt1+pt2)*T0p(iir,iit,iip) &
                                    !& + c(iir,iit,iip)*(pp1+pp2)*T0p(iir,iit,iip) - f(iir,iit,iip)*(pp1+pp2)*T0t(iir,iit,iip) &
                                    !& + a(iir,iit,iip)*T0r(iir,iit,iip)**2 &
                                    !& + b(iir,iit,iip)*T0t(iir,iit,iip)**2 + c(iir,iit,iip)*T0p(iir,iit,iip)**2 &
                                    !& - 2*f(iir,iit,iip)*T0t(iir,iit,iip)*T0p(iir,iit,iip) )
                                    
                                    Htau = sqrt( &
                                    &   a(iir,iit,iip)*(T0r(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pr1+pr2)/2)**2 &
                                    & + b(iir,iit,iip)*(T0t(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pt1+pt2)/2)**2 &
                                    & + c(iir,iit,iip)*(T0p(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pp1+pp2)/2)**2 &
                                    &-2*f(iir,iit,iip)*(T0t(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pt1+pt2)/2) &
                                    &                 *(T0p(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pp1+pp2)/2) )
                                    

                                    ! 更新 timetable
                                    tpT=coe*(fun(iir,iit,iip)-Htau)  &
                                    & +coe*(sigr*(pr2-pr1)/2+sigt*(pt2-pt1)/2+sigp*(pp2-pp1)/2)+tau(iir,iit,iip);

                                    !write(*,*) fun(iir,iit,iip),Htau

                                    tau(iir,iit,iip) = tpT

                                    
                                    mpi_ptn=mpi_ptn+1
                                    mpi_ptijk(:,mpi_ptn) = (/iir,iit,iip/)
                                    mpi_ptv(mpi_ptn) = tpT
                                    my_npt = my_npt + 1
                                end if                              
                            end do
                        end do
                        
                        !if (1>2) then
                        ! 数据通信 步骤1, 其他线程传输给主线程
                        if (myid .eq. 0) then
                            ! 负责接受数据
                            do iproc=1,nproc-1
                                call mpi_recv(int_temp,1,mpi_integer,iproc,tag,mpi_comm_world,istat,ierr)
                                if (int_temp .ne. 0) then
                                    call mpi_recv(mpi_ptijk(1,mpi_ptn+1),int_temp*3,mpi_integer,iproc,tag+1, &
                                                & mpi_comm_world,istat,ierr)
                                    call mpi_recv(mpi_ptv(mpi_ptn+1),int_temp,mpi_double_precision,iproc,tag+2, &
                                                & mpi_comm_world,istat,ierr)
                                    ! 接收完数据，赋值
                                    do ipt=mpi_ptn+1,mpi_ptn+int_temp
                                        tau(mpi_ptijk(1,ipt),mpi_ptijk(2,ipt),mpi_ptijk(3,ipt)) = mpi_ptv(ipt)
                                    end do
                                    mpi_ptn = mpi_ptn + int_temp
                                end if
                            end do
                        else
                            ! 负责发送数据
                            call mpi_send(mpi_ptn,1,mpi_integer,0,tag,mpi_comm_world,ierr)
                            !print *, mpi_ptn, mpi_ptijk(:,1:mpi_ptn), mpi_ptv(1:mpi_ptn)
                            if (mpi_ptn .ne. 0) then
                                call mpi_send(mpi_ptijk,3*mpi_ptn,mpi_integer,0,tag+1,mpi_comm_world,ierr)
                                call mpi_send(mpi_ptv,mpi_ptn,mpi_double_precision,0,tag+2,mpi_comm_world,ierr)    
                            end if                       
                        end if
                        
                        ! 数据通信 步骤2, 主线程将更新后的level数据广播到其他线程

                        call mpi_bcast(mpi_ptn,1,mpi_integer,0,mpi_comm_world,ierr)
                        call mpi_bcast(mpi_ptijk,3*mpi_ptn,mpi_integer,0,mpi_comm_world,ierr)
                        call mpi_bcast(mpi_ptv,mpi_ptn,mpi_double_precision,0,mpi_comm_world,ierr)    
                        
                        ! 数据广播完成，开始赋值
                        if (myid .ne. 0) then
                            do ipt=1,mpi_ptn
                                tau(mpi_ptijk(1,ipt),mpi_ptijk(2,ipt),mpi_ptijk(3,ipt)) = mpi_ptv(ipt)
                            end do
                        end if

                        call mpi_barrier(mpi_comm_world,ierr)
                        !end if
                    end do    

                    ! print *, my_npt, Nptn
                    
                    ! 处理边界

                    do iit=1,nt
                        do iip=1,np
                            if (tau(3,iit,iip)>0) then
                                tau(1,iit,iip) = max(2*tau(2,iit,iip)-tau(3,iit,iip),tau(3,iit,iip))
                            else
                                tau(1,iit,iip) = min(2*tau(2,iit,iip)-tau(3,iit,iip),tau(3,iit,iip))
                            end if
                            if (tau(nr-2,iit,iip)>0) then
                                tau(nr,iit,iip) = max(2*tau(nr-1,iit,iip)-tau(nr-2,iit,iip),tau(nr-2,iit,iip))
                            else
                                tau(nr,iit,iip) = min(2*tau(nr-1,iit,iip)-tau(nr-2,iit,iip),tau(nr-2,iit,iip))
                            end if
                        end do
                    end do
                    do iir=1,nr
                        do iip=1,np
                            if (tau(iir,3,iip)>0) then
                                tau(iir,1,iip) = max(2*tau(iir,2,iip)-tau(iir,3,iip),tau(iir,3,iip))
                            else
                                tau(iir,1,iip) = min(2*tau(iir,2,iip)-tau(iir,3,iip),tau(iir,3,iip))
                            end if
                            if (tau(iir,nt-2,iip)>0) then
                                tau(iir,nt,iip) = max(2*tau(iir,nt-1,iip)-tau(iir,nt-2,iip),tau(iir,nt-2,iip))
                            else
                                tau(iir,nt,iip) = min(2*tau(iir,nt-1,iip)-tau(iir,nt-2,iip),tau(iir,nt-2,iip))
                            end if
                        end do
                    end do
                    do iir=1,nr
                        do iit=1,nt
                            if (tau(iir,iit,3)>0) then
                                tau(iir,iit,1) = max(2*tau(iir,iit,2)-tau(iir,iit,3),tau(iir,iit,3))
                            else 
                                tau(iir,iit,1) = min(2*tau(iir,iit,2)-tau(iir,iit,3),tau(iir,iit,3))
                            end if
                            if (tau(iir,iit,np-2)>0) then
                                tau(iir,iit,np) = max(2*tau(iir,iit,np-1)-tau(iir,iit,np-2),tau(iir,iit,np-2))
                            else
                                tau(iir,iit,np) = min(2*tau(iir,iit,np-1)-tau(iir,iit,np-2),tau(iir,iit,np-2))
                            end if
                        end do
                    end do   


                end do
            end do
        end do

        

        ! 统计误差，判断迭代终止条件
        do iir=2,nr-1
            do iit=2,nt-1
                do iip=2,np-1
                    L1_dif=L1_dif+abs(tau(iir,iit,iip)-tau_old(iir,iit,iip))*T0v(iir,iit,iip)
                    Linf_dif=max(Linf_dif,abs(tau(iir,iit,iip)-tau_old(iir,iit,iip))*T0v(iir,iit,iip))
                end do               
            end do
        end do
        
        do iir=3,nr-2
            do iit=3,nt-2
                do iip=3,np-2
                    L1_err=L1_err+abs(u(iir,iit,iip)-T0v(iir,iit,iip)*tau(iir,iit,iip))
                    Linf_err=max(Linf_err,abs(T0v(iir,iit,iip)*tau(iir,iit,iip)-u(iir,iit,iip)))   
                end do            
            end do
        end do
        L1_err=L1_err/((nr-4)*(nt-4)*(np-4))
        L1_dif=L1_dif/((nr-2)*(nt-2)*(np-2))

        if (myid .eq. 0) then 
            ! ################  iteration information  #################
            if (abs(L1_dif)<tol) then
                write(*,*) 'id',myid,'iter ',iter,', T is steadt'
                isexit = .true.
            else
                write(*,*) 'id',myid,'iter ',iter,', T is changing, L1 dif = ', L1_dif,'L inf dif = ', Linf_dif
            end if

            if (iter==MaxIter) then    
                write(*,*) 'id','iter ',iter,', max iteration steps'
            end if

            ! ################  solver accuracy  #################
            !write(*,'(a,f15.7,a,f15.7)') 'L_1(T-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err
        end if

        call mpi_bcast(isexit,1,mpi_logical,0,mpi_comm_world,ierr)

        if (isexit) then
            exit
        end if
    end do

    T=T0v*tau
    deallocate(proc_irange)
end subroutine



subroutine FSM_O1_PS_sphe_3d_mul_mpi(rr,tt,pp,nr,nt,np,spha,sphb,sphc,sphf,T,fun,r0,t0,p0,u)
    
    use mpi
    ! a -d -e
    ! -d b -f
    ! -e -f c
    integer nr,nt,np
    double precision :: dr,dt,dp,rr(nr),tt(nt),pp(np),a(nr,nt,np),b(nr,nt,np),c(nr,nt,np),f(nr,nt,np)
    double precision :: spha(nr,nt,np),sphb(nr,nt,np),sphc(nr,nt,np),sphf(nr,nt,np)
    double precision :: fun(nr,nt,np),T(nr,nt,np),ischange(nr,nt,np),u(nr,nt,np),xi(nr,nt,np),eta(nr,nt,np)
    double precision :: r1,r2,r3,a0,b0,c0,d0,f0,T0v(nr,nt,np),T0r(nr,nt,np),T0t(nr,nt,np),T0p(nr,nt,np)
    double precision :: tau(nr,nt,np),tau_old(nr,nt,np),r0,t0,p0
    integer :: idr0,idt0,idp0
    double precision :: sigr,sigt,sigp,coe,pr1,pr2,wr1,wr2,pt1,pt2,wt1,wt2,pp1,pp2,wp1,wp2,tpT
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err
    double precision :: t1,p1
    double precision,parameter :: tol = (10.0)**(-3),eps=10.0**(-12)
    integer,parameter :: MaxIter=200
    integer iter,rdirec,tdirec,pdirec,iir,iit,iip
    double precision :: x0,y0,z0,x,y,z,xst,yst,zst,e11,e12,e13,e21,e22,e23,e31,e32,e33
    double precision :: xstr,xstt,xstp,ystr,ystt,ystp,zstr,zstt,zstp
    integer :: ileft,iright,jleft,jright,ii,jj,kk
    logical :: isexit

    integer :: ierr,myid,nproc,tag,istat(mpi_status_size),iproc,ptid(nr,np,nt)
    integer :: mpi_ptijk(3,np*nt),mpi_ptn,ipt,int_temp,my_npt
    double precision :: mpi_ptv(np*nt),dp_temp

    integer :: ptr1,ptr2,iNumber,total_ptn,Nptn,sNptn,tpNptn
    double precision :: ave_ni
    integer,allocatable :: proc_irange(:,:,:)


    isexit = .false.
    tag=99
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    
    ! ------------------------ 构造网格 ------------------------ 
    dr=rr(2)-rr(1); dt=tt(2)-tt(1); dp=pp(2)-pp(1)

    ! ------------------------ 构造矩阵 -------------------------  
    !    a -d -e
    !   -d  b -f
    !   -e -f  c

    ! ------------------------ 构造 T0 ------------------------ 
    
    ! 震源处参数离散化
    idr0=floor((r0-rr(1))/dr+1); idt0=floor((t0-tt(1))/dt+1); idp0=floor((p0-pp(1))/dp+1); 
    r1 = min(1.0,(r0-rr(idr0))/dr); r2 = min(1.0,(t0-tt(idt0))/dt); r3 = min(1.0,(p0-pp(idp0))/dp); 

    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                a(iir,iit,iip) = spha(iir,iit,iip)
                b(iir,iit,iip) = sphb(iir,iit,iip)/(rr(iir)**2)
                c(iir,iit,iip) = sphc(iir,iit,iip)/(rr(iir)**2*cos(tt(iit))**2)
                f(iir,iit,iip) = sphf(iir,iit,iip)/(rr(iir)**2*cos(tt(iit)))
            end do
        end do
    end do
    
    a0=(1-r1)*(1-r2)*(1-r3)*a(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*a(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*a(idr0,idt0+1,idp0)+(1-r1)*r2*r3*a(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*a(idr0+1,idt0,idp0)+r1*(1-r2)*r3*a(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*a(idr0+1,idt0+1,idp0)+r1*r2*r3*a(idr0+1,idt0+1,idp0+1)

    b0=(1-r1)*(1-r2)*(1-r3)*b(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*b(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*b(idr0,idt0+1,idp0)+(1-r1)*r2*r3*b(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*b(idr0+1,idt0,idp0)+r1*(1-r2)*r3*b(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*b(idr0+1,idt0+1,idp0)+r1*r2*r3*b(idr0+1,idt0+1,idp0+1)

    c0=(1-r1)*(1-r2)*(1-r3)*c(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*c(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*c(idr0,idt0+1,idp0)+(1-r1)*r2*r3*c(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*c(idr0+1,idt0,idp0)+r1*(1-r2)*r3*c(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*c(idr0+1,idt0+1,idp0)+r1*r2*r3*c(idr0+1,idt0+1,idp0+1)

    f0=(1-r1)*(1-r2)*(1-r3)*f(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*f(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*f(idr0,idt0+1,idp0)+(1-r1)*r2*r3*f(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*f(idr0+1,idt0,idp0)+r1*(1-r2)*r3*f(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*f(idr0+1,idt0+1,idp0)+r1*r2*r3*f(idr0+1,idt0+1,idp0+1)
    !d0=-d0

    fun0=(1-r1)*(1-r2)*(1-r3)*fun(idr0,idt0,idp0)+(1-r1)*(1-r2)*r3*fun(idr0,idt0,idp0+1) &
    & +(1-r1)*r2*(1-r3)*fun(idr0,idt0+1,idp0)+(1-r1)*r2*r3*fun(idr0,idt0+1,idp0+1) &
    & +r1*(1-r2)*(1-r3)*fun(idr0+1,idt0,idp0)+r1*(1-r2)*r3*fun(idr0+1,idt0,idp0+1) &
    & +r1*r2*(1-r3)*fun(idr0+1,idt0+1,idp0)+r1*r2*r3*fun(idr0+1,idt0+1,idp0+1)

    

    ! 构造T0
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                r1 = rr(iir)
                t1 = tt(iit)
                p1 = pp(iip)
       
                T0v(iir,iit,iip) = fun0*sqrt( 1.0/a0*(r1-r0)**2 + c0/(c0*b0-f0**2)*(t1-t0)**2 &
                                          & + b0/(c0*b0-f0**2)*(p1-p0)**2 + 2*f0/(c0*b0-f0**2)*(t1-t0)*(p1-p0) ) 
                ischange(iir,iit,iip)=1
                if ( T0v(iir,iit,iip)==0 ) then
                    T0r(iir,iit,iip) = 0
                    T0t(iir,iit,iip) = 0
                    T0p(iir,iit,iip) = 0
                else
                    T0r(iir,iit,iip) = fun0**2*(1.0/a0*(r1-r0))/T0v(iir,iit,iip)
                    T0t(iir,iit,iip) = fun0**2*(c0/(c0*b0-f0**2)*(t1-t0)+f0/(c0*b0-f0**2)*(p1-p0))/T0v(iir,iit,iip)
                    T0p(iir,iit,iip) = fun0**2*(b0/(c0*b0-f0**2)*(p1-p0)+f0/(c0*b0-f0**2)*(t1-t0))/T0v(iir,iit,iip)    
                end if


                if ( abs((rr(iir)-r0)/dr)<=2 .and. abs((tt(iit)-t0)/dt)<=2 .and. abs((pp(iip)-p0)/dp)<= 2) then
                    tau(iir,iit,iip) = 1  !震源周围几个点，直接认为是常速度结构，给出解析解，即和T0相等
                    !tau(iir,iit,iip) = u(iir,iit,iip) - T0v(iir,iit,iip)
                    ischange(iir,iit,iip)=0
                    if (iir==1 .or. iir==nr .or. iit==1 .or. iit==nt .or. iip==1 .or. iip==np) then
                        write(*,*) 'source on the boundary, mesh error'
                        print *, rr(iir),t0*180/3.1415927,iit,p0*180/3.1415927,iip
                    !    pause
                    end if
                    !write(*,*) iir-idr0,iit-idt0,iip-idp0,u(iir,iit,iip),T0v(iir,iit,iip),u(iir,iit,iip)-T0v(iir,iit,iip)
                else
                    tau(iir,iit,iip) = 1
                    ischange(iir,iit,iip)=1
                end if
            end do
        end do
    end do


    L1_err=0; Linf_err=0
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                L1_err=L1_err+abs(u(iir,iit,iip)-T0v(iir,iit,iip))
                Linf_err=max(Linf_err,abs(T0v(iir,iit,iip)-u(iir,iit,iip)))   
                if(isnan(L1_err)) then
                    print *, u(iir,iit,iip),T0v(iir,iit,iip)
                    pause
                end if
            end do            
        end do
    end do
    L1_err = L1_err/(nr*np*nt)
    !if (myid .eq. 0) then
    !    write(*,*) 'T0 L1 error is ',L1_err,', T0 Linf error is', Linf_err
    !end if

    allocate(proc_irange(2,nproc,nr+nt+np-3-6+1))

    ! 给每个线程分配计算点
    Nptn=0
    sNptn=0

    do level = 6,nr+nt+np-3
        ileft = max(2,level-np-nt+2)
        iright = min(level-4,nr-1)

        iNumber = iright - ileft + 1

        total_ptn = 0
        do ip = ileft,iright
            total_ptn = total_ptn+(min(level-ip-2,nt-1)-max(2,level-ip-np+1)+1)
        end do
        Nptn = Nptn + total_ptn  ! 总点数
        ave_ni = total_ptn*1.0/nproc

        if (iNumber <= nproc) then
            ! 即i的数量小于等于线程数，每个线程分配一个i就好
            tpNptn = 0
            MaxtpNptn = 0
            do ip = ileft,iright
                proc_irange(:,ip-ileft+1,level-5)=(/ip,ip/)
                tpNptn = (min(level-ip-2,nt-1)-max(2,level-ip-np+1)+1)
                MaxtpNptn = max(MaxtpNptn,tpNptn) 
            end do
            do ip = iright+1,ileft+nproc-1
                proc_irange(:,ip-ileft+1,level-5)=(/-1,-2/)
            end do
            sNptn = sNptn + MaxtpNptn
        else
            !if (myid .eq. 0) then
            !    print *, ileft,iright
            !end if
            !call sleep(2)
            
            ! i的数量大于线程数，平均分配，接近 total_ptn/nproc
            int_temp = 0
            tpNptn = 0
            MaxtpNptn = 0
            ptr1 = ileft
            iproc = 1
            do ptr2 = ileft,iright             
                int_temp = int_temp + (min(level-ptr2-2,nt-1)-max(2,level-ptr2-np+1)+1)
                tpNptn = tpNptn + (min(level-ptr2-2,nt-1)-max(2,level-ptr2-np+1)+1)
                !if (myid .eq. 0) then
                !    print *, ptr2,(min(level-ptr2-2,nt-1)-max(2,level-ptr2-np+1)+1),int_temp,ave_ni
                !end if
                !call sleep(2)
                if ((int_temp>=ave_ni*iproc) .or. (ptr2 .eq.iright)) then
                    !if (iproc>nproc) then
                    !    print *, iproc,nproc
                    !end if
                    proc_irange(:,iproc,level-5) = (/ptr1,ptr2/)
                    ptr1 = ptr2+1
                    iproc = iproc +1
                    MaxtpNptn = max(MaxtpNptn,tpNptn) 
                    tpNptn = 0
                end if
            end do   
            sNptn = sNptn + MaxtpNptn         
            !call sleep(2)
        end if
        !if (myid .eq. 0) then
        !    print *, proc_irange(1,:,level-5)
        !    print *, proc_irange(2,:,level-5)
        !    print *, ' '
        !end if
    end do

    ! print *, Nptn*1.0/sNptn


    ! 正式迭代，更新tau
    do iter =1,MaxIter
        tau_old = tau
        L1_dif=0; Linf_dif=0;L1_err=0;Linf_err=0;
        
        do rdirec = -1,1,2
            do tdirec = -1,1,2
                do pdirec = -1,1,2
                    my_npt = 0
                    do level = 6,nr+nt+np-3
                        ! 2<= ir <= nr-1;   2<= it <= nr-1;   2<=ip <=np-1

                        ileft = max(2,level-np-nt+2)
                        iright = min(level-4,nr-1)

                        mpi_ptn=0
                        do ii = proc_irange(1,myid+1,level-5),proc_irange(2,myid+1,level-5)
                        
                            jleft = max(2,level-ii-np+1)
                            jright = min(level-ii-2,nt-1)
                            do jj = jleft,jright

                                kk = level-ii-jj
                                
                                if(rdirec<0) then
                                    iir=ii
                                else
                                    iir=nr+1-ii
                                end if
                                if(tdirec<0) then
                                    iit=jj
                                else
                                    iit=nt+1-jj
                                end if
                                if(pdirec<0) then
                                    iip=kk
                                else
                                    iip=np+1-kk
                                end if
            

                                if(ischange(iir,iit,iip)==1) then
                                    

                                    !sigr=2*sqrt(a(iir,iit,iip)); sigt=2*sqrt(b(iir,iit,iip)); sigp=2*sqrt(c(iir,iit,iip))
                                    !coe=1.0/((sigr/dr)+(sigt/dt)+(sigp/dp))

                                    sigr = 2*sqrt(a(iir,iit,iip))*T0v(iir,iit,iip)
                                    sigt = 2*sqrt(b(iir,iit,iip))*T0v(iir,iit,iip)
                                    sigp = 2*sqrt(c(iir,iit,iip))*T0v(iir,iit,iip)
                                    coe=1.0/((sigr/dr)+(sigt/dt)+(sigp/dp))

                                    ! 构造单方向梯度 1 阶格式
                                    pr1=(tau(iir,iit,iip)-tau(iir-1,iit,iip))/dr;
                                    pr2=(tau(iir+1,iit,iip)-tau(iir,iit,iip))/dr;          
                                    
                                    pt1=(tau(iir,iit,iip)-tau(iir,iit-1,iip))/dt; 
                                    pt2=(tau(iir,iit+1,iip)-tau(iir,iit,iip))/dt;  
                                    
                                    pp1=(tau(iir,iit,iip)-tau(iir,iit,iip-1))/dp;
                                    pp2=(tau(iir,iit,iip+1)-tau(iir,iit,iip))/dp; 

                                    !计算 LF Hamiltonian 

                                    ! additive
                                    !Htau=sqrt( a(iir,iit,iip)*((pr1+pr2)/2)**2 + b(iir,iit,iip)*((pt1+pt2)/2)**2 &
                                    !& + c(iir,iit,iip)*((pp1+pp2)/2)**2 -2*f(iir,iit,iip)*(pt1+pt2)/2*(pp1+pp2)/2 &
                                    !& + a(iir,iit,iip)*(pr1+pr2)*T0r(iir,iit,iip) &
                                    !& + b(iir,iit,iip)*(pt1+pt2)*T0t(iir,iit,iip) - f(iir,iit,iip)*(pt1+pt2)*T0p(iir,iit,iip) &
                                    !& + c(iir,iit,iip)*(pp1+pp2)*T0p(iir,iit,iip) - f(iir,iit,iip)*(pp1+pp2)*T0t(iir,iit,iip) &
                                    !& + a(iir,iit,iip)*T0r(iir,iit,iip)**2 &
                                    !& + b(iir,iit,iip)*T0t(iir,iit,iip)**2 + c(iir,iit,iip)*T0p(iir,iit,iip)**2 &
                                    !& - 2*f(iir,iit,iip)*T0t(iir,iit,iip)*T0p(iir,iit,iip) )   
                                    
                                    ! multiplicative
                                    Htau = sqrt( &
                                    &   a(iir,iit,iip)*(T0r(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pr1+pr2)/2)**2 &
                                    & + b(iir,iit,iip)*(T0t(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pt1+pt2)/2)**2 &
                                    & + c(iir,iit,iip)*(T0p(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pp1+pp2)/2)**2 &
                                    &-2*f(iir,iit,iip)*(T0t(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pt1+pt2)/2) &
                                    &                 *(T0p(iir,iit,iip)*tau(iir,iit,iip)+T0v(iir,iit,iip)*(pp1+pp2)/2) )
                                    

                                    ! 更新 timetable
                                    tpT=coe*(fun(iir,iit,iip)-Htau)  &
                                    & +coe*(sigr*(pr2-pr1)/2+sigt*(pt2-pt1)/2+sigp*(pp2-pp1)/2)+tau(iir,iit,iip);

                                    !write(*,*) fun(iir,iit,iip),Htau

                                    tau(iir,iit,iip) = tpT

                                    
                                    mpi_ptn=mpi_ptn+1
                                    mpi_ptijk(:,mpi_ptn) = (/iir,iit,iip/)
                                    mpi_ptv(mpi_ptn) = tpT
                                    my_npt = my_npt + 1
                                end if                              
                            end do
                        end do
                        
                        !if (1>2) then
                        ! 数据通信 步骤1, 其他线程传输给主线程
                        if (myid .eq. 0) then
                            ! 负责接受数据
                            do iproc=1,nproc-1
                                call mpi_recv(int_temp,1,mpi_integer,iproc,tag,mpi_comm_world,istat,ierr)
                                if (int_temp .ne. 0) then
                                    call mpi_recv(mpi_ptijk(1,mpi_ptn+1),int_temp*3,mpi_integer,iproc,tag+1, &
                                                & mpi_comm_world,istat,ierr)
                                    call mpi_recv(mpi_ptv(mpi_ptn+1),int_temp,mpi_double_precision,iproc,tag+2, &
                                                & mpi_comm_world,istat,ierr)
                                    ! 接收完数据，赋值
                                    do ipt=mpi_ptn+1,mpi_ptn+int_temp
                                        tau(mpi_ptijk(1,ipt),mpi_ptijk(2,ipt),mpi_ptijk(3,ipt)) = mpi_ptv(ipt)
                                    end do
                                    mpi_ptn = mpi_ptn + int_temp
                                end if
                            end do
                        else
                            ! 负责发送数据
                            call mpi_send(mpi_ptn,1,mpi_integer,0,tag,mpi_comm_world,ierr)
                            !print *, mpi_ptn, mpi_ptijk(:,1:mpi_ptn), mpi_ptv(1:mpi_ptn)
                            if (mpi_ptn .ne. 0) then
                                call mpi_send(mpi_ptijk,3*mpi_ptn,mpi_integer,0,tag+1,mpi_comm_world,ierr)
                                call mpi_send(mpi_ptv,mpi_ptn,mpi_double_precision,0,tag+2,mpi_comm_world,ierr)    
                            end if                       
                        end if
                        
                        ! 数据通信 步骤2, 主线程将更新后的level数据广播到其他线程

                        call mpi_bcast(mpi_ptn,1,mpi_integer,0,mpi_comm_world,ierr)
                        call mpi_bcast(mpi_ptijk,3*mpi_ptn,mpi_integer,0,mpi_comm_world,ierr)
                        call mpi_bcast(mpi_ptv,mpi_ptn,mpi_double_precision,0,mpi_comm_world,ierr)    
                        
                        ! 数据广播完成，开始赋值
                        if (myid .ne. 0) then
                            do ipt=1,mpi_ptn
                                tau(mpi_ptijk(1,ipt),mpi_ptijk(2,ipt),mpi_ptijk(3,ipt)) = mpi_ptv(ipt)
                            end do
                        end if

                        call mpi_barrier(mpi_comm_world,ierr)
                        !end if
                    end do    

                    ! print *, my_npt, Nptn
                    
                    ! 处理边界

                    do iit=1,nt
                        do iip=1,np
                            tau(1,iit,iip) = max(2*tau(2,iit,iip)-tau(3,iit,iip),tau(3,iit,iip))
                            tau(nr,iit,iip) = max(2*tau(nr-1,iit,iip)-tau(nr-2,iit,iip),tau(nr-2,iit,iip))
                        end do
                    end do
                    do iir=1,nr
                        do iip=1,np
                            tau(iir,1,iip) = max(2*tau(iir,2,iip)-tau(iir,3,iip),tau(iir,3,iip))                       
                            tau(iir,nt,iip) = max(2*tau(iir,nt-1,iip)-tau(iir,nt-2,iip),tau(iir,nt-2,iip))
                        end do
                    end do
                    do iir=1,nr
                        do iit=1,nt
                            tau(iir,iit,1) = max(2*tau(iir,iit,2)-tau(iir,iit,3),tau(iir,iit,3))                          
                            tau(iir,iit,np) = max(2*tau(iir,iit,np-1)-tau(iir,iit,np-2),tau(iir,iit,np-2))
                        end do
                    end do   

                end do
            end do
        end do

        

        ! 统计误差，判断迭代终止条件
        do iir=2,nr-1
            do iit=2,nt-1
                do iip=2,np-1
                    L1_dif=L1_dif+abs(tau(iir,iit,iip)-tau_old(iir,iit,iip))*T0v(iir,iit,iip)
                    Linf_dif=max(Linf_dif,abs(tau(iir,iit,iip)-tau_old(iir,iit,iip))*T0v(iir,iit,iip))
                end do               
            end do
        end do
        
        do iir=3,nr-2
            do iit=3,nt-2
                do iip=3,np-2
                    L1_err=L1_err+abs(u(iir,iit,iip)-T0v(iir,iit,iip)*tau(iir,iit,iip))
                    Linf_err=max(Linf_err,abs(T0v(iir,iit,iip)*tau(iir,iit,iip)-u(iir,iit,iip)))   
                end do            
            end do
        end do
        L1_err=L1_err/((nr-4)*(nt-4)*(np-4))
        L1_dif=L1_dif/((nr-2)*(nt-2)*(np-2))

        if (myid .eq. 0) then 
            ! ################  iteration information  #################
            if (abs(L1_dif)<tol) then
                write(*,*) 'id',myid,'iter ',iter,', T is steadt'
                isexit = .true.
            else
                write(*,*) 'id',myid,'iter ',iter,', T is changing, L1 dif = ', L1_dif,'L inf dif = ', Linf_dif
            end if

            if (iter==MaxIter) then    
                write(*,*) 'id','iter ',iter,', max iteration steps'
            end if

            ! ################  solver accuracy  #################
            !write(*,'(a,f15.7,a,f15.7)') 'L_1(T-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err
        end if

        call mpi_bcast(isexit,1,mpi_logical,0,mpi_comm_world,ierr)

        if (isexit) then
            exit
        end if
    end do

    T=T0v*tau
    
end subroutine

subroutine FSM_O1_sphe_3d_mpi(rr,tt,pp,nr,nt,np,spha,sphb,sphc,sphf,T,fun,r0,t0,p0,u)
    use mpi
    ! a -d -e
    ! -d b -f
    ! -e -f c
    integer nr,nt,np
    double precision :: dr,dt,dp,rr(nr),tt(nt),pp(np),a(nr,nt,np),b(nr,nt,np),c(nr,nt,np),f(nr,nt,np)
    double precision :: spha(nr,nt,np),sphb(nr,nt,np),sphc(nr,nt,np),sphf(nr,nt,np)
    double precision :: fun(nr,nt,np),T(nr,nt,np),ischange(nr,nt,np),u(nr,nt,np),xi(nr,nt,np),eta(nr,nt,np)
    double precision :: r1,r2,r3,a0,b0,c0,d0,f0,T0v(nr,nt,np),T0r(nr,nt,np),T0t(nr,nt,np),T0p(nr,nt,np)
    double precision :: tau(nr,nt,np),tau_old(nr,nt,np),r0,t0,p0
    integer :: idr0,idt0,idp0
    double precision :: sigr,sigt,sigp,coe,pr1,pr2,wr1,wr2,pt1,pt2,wt1,wt2,pp1,pp2,wp1,wp2,tpT
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err
    double precision :: t1,p1
    double precision,parameter :: tol = (10.0)**(-3),eps=10.0**(-12)
    integer,parameter :: MaxIter=1000
    integer iter,rdirec,tdirec,pdirec,iir,iit,iip
    double precision :: x0,y0,z0,x,y,z,xst,yst,zst,e11,e12,e13,e21,e22,e23,e31,e32,e33
    double precision :: xstr,xstt,xstp,ystr,ystt,ystp,zstr,zstt,zstp
    integer :: ileft,iright,jleft,jright,ii,jj,kk
    logical :: isexit

    integer :: ierr,myid,nproc,tag,istat(mpi_status_size),iproc,ptid(nr,np,nt)
    integer :: mpi_ptijk(3,np*nt),mpi_ptn,ipt,int_temp,my_npt
    double precision :: mpi_ptv(np*nt),dp_temp

    integer :: ptr1,ptr2,iNumber,total_ptn,Nptn,sNptn,tpNptn
    double precision :: ave_ni
    integer,allocatable :: proc_irange(:,:,:)


    isexit = .false.
    tag=99
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    call mpi_comm_size(mpi_comm_world,nproc,ierr)

    
    ! ------------------------ 构造网格 build mesh ------------------------ 
    dr=rr(2)-rr(1); dt=tt(2)-tt(1); dp=pp(2)-pp(1)

    ! ------------------------ 构造矩阵 -------------------------  
    !    a -d -e
    !   -d  b -f
    !   -e -f  c
    
    ! 震源处参数离散化 source parameterization
    idr0=floor((r0-rr(1))/dr+1); idt0=floor((t0-tt(1))/dt+1); idp0=floor((p0-pp(1))/dp+1); 
    r1 = min(1.0,(r0-rr(idr0))/dr); r2 = min(1.0,(t0-tt(idt0))/dt); r3 = min(1.0,(p0-pp(idp0))/dp); 

    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                a(iir,iit,iip) = spha(iir,iit,iip)
                b(iir,iit,iip) = sphb(iir,iit,iip)/(rr(iir)**2)
                c(iir,iit,iip) = sphc(iir,iit,iip)/(rr(iir)**2*cos(tt(iit))**2)
                f(iir,iit,iip) = sphf(iir,iit,iip)/(rr(iir)**2*cos(tt(iit)))
            end do
        end do
    end do

    ! initla value of tau
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                if ( ((rr(iir)-r0)/dr)**2 + ((tt(iit)-t0)/dt)**2 + ((pp(iip)-p0)/dp)**2 < 2**2) then
                    tau(iir,iit,iip) = 0    ! source condition
                    ischange(iir,iit,iip)=0
                    if (iir==1 .or. iir==nr .or. iit==1 .or. iit==nt .or. iip==1 .or. iip==np) then
                        write(*,*) 'source on the boundary, mesh error'
                        print *, rr(iir),r0,tt(iit),t0,pp(iip),p0
                        pause
                    end if
                else
                    tau(iir,iit,iip) = 300
                    ischange(iir,iit,iip)=1
                end if
            end do
        end do
    end do


    allocate(proc_irange(2,nproc,nr+nt+np-3-6+1))

    ! we have origin 3 loops for i,j,k:   
    !   do i=2,nr-1
    !       do j=2,nt-1
    !           do k=2,np-1
    !               calculate point x(i,j,k) 

    ! we can also loop level,i,j. Here level =i+j+k, ranging from 2+2+2 to nx-1+ny-1+nz-1, that is, [6, nx+ny+nz-3]   
    !   do level = 6, nr+nt+np-3            
    !       do i=ileft,iright               ileft = max(2,level-np-nt+2); iright = min(level-4,nr-1)
    !           do j=jleft,jright           jleft = max(2,level-ii-np+1); jright = min(level-ii-2,nt-1)
    !               calculate point x(i,j,k) = x(i,j,level-i-j)               

    ! we know that we can parallelize the calculation at each points in the same level.
    ! Here we choose to parallelize the loop i.
    ! e.g., for a fixed level , in processor 'n', we only calculate points:
    !       do i = proc_irange(1,n,level-5), proc_irange(2,n,level-5).   
    !           do j =jleft,jright
    !               calculate point x(i,j,k) = x(i,j,level-i-j)

    ! Here we define proc_irange(m,n,p). m = 1 or 2 means the left or right boundary of index i.
    !                     n is the id of the processor,  p represent the level-5, ranging from 1 to nx+ny+nz-8 (level = 6,7,8, ... nx+ny+nz-3)

    ! what we need is to calculate proc_irange(m,n,p) so that for each processor, to number of points are similar (load balancing)     

    ! --------------- load balancing --------------
    Nptn=0
    sNptn=0

    do level = 6,nr+nt+np-3
        ! ----- for each level -----
        ileft = max(2,level-np-nt+2)        ! we know the range of index i = ileft, iright
        iright = min(level-4,nr-1)          

        iNumber = iright - ileft + 1        ! the number of i 

        total_ptn = 0       ! the total number of points for the fixed level 
        do ip = ileft,iright
            ! for each index i, we sum the number of points in loop j=jleft,jright, where jleft = max(2,level-ii-np+1); jright = min(level-ii-2,nt-1)
            total_ptn = total_ptn+(min(level-ip-2,nt-1)-max(2,level-ip-np+1)+1) 
        end do
        Nptn = Nptn + total_ptn  ! the total number of points
        ave_ni = total_ptn*1.0/nproc

        if (iNumber <= nproc) then
            ! we parallelize loop i. if number of index i is smaller than the number of processor, each processor calculate one index i
            tpNptn = 0
            MaxtpNptn = 0
            do ip = ileft,iright
                proc_irange(:,ip-ileft+1,level-5)=(/ip,ip/)     ! processor 'ip-ileft+1' calculate index 'ip'   
                tpNptn = (min(level-ip-2,nt-1)-max(2,level-ip-np+1)+1)
                MaxtpNptn = max(MaxtpNptn,tpNptn) 
            end do
            do ip = iright+1,ileft+nproc-1      !  no index i for the rest processors, just giving invalid loof range.
                proc_irange(:,ip-ileft+1,level-5)=(/-1,-2/)
            end do
            sNptn = sNptn + MaxtpNptn
        else            
            ! i的数量大于线程数，平均分配，接近 total_ptn/nproc
            ! if ithe number of index i is greater than the number of processor, we try to make that the number of points 
            ! in each processor is close to the average, that is total_ptn/nproc ( total number / number of processors)
            int_temp = 0
            tpNptn = 0
            MaxtpNptn = 0
            ptr1 = ileft
            iproc = 1
            do ptr2 = ileft,iright             
                int_temp = int_temp + (min(level-ptr2-2,nt-1)-max(2,level-ptr2-np+1)+1)
                tpNptn = tpNptn + (min(level-ptr2-2,nt-1)-max(2,level-ptr2-np+1)+1)

                if ((int_temp>=ave_ni*iproc) .or. (ptr2 .eq.iright)) then

                    proc_irange(:,iproc,level-5) = (/ptr1,ptr2/)
                    ptr1 = ptr2+1
                    iproc = iproc +1
                    MaxtpNptn = max(MaxtpNptn,tpNptn) 
                    tpNptn = 0
                end if
            end do   
            sNptn = sNptn + MaxtpNptn         
        end if
    end do

    ! print *, Nptn*1.0/sNptn


    ! 正式迭代，更新tau
    do iter =1,MaxIter
        tau_old = tau
        L1_dif=0; Linf_dif=0;L1_err=0;Linf_err=0;
        
        do rdirec = -1,1,2
            do tdirec = -1,1,2
                do pdirec = -1,1,2
                    my_npt = 0
                    do level = 6,nr+nt+np-3
                        ! 2<= ir <= nr-1;   2<= it <= nr-1;   2<=ip <=np-1

                        ileft = max(2,level-np-nt+2)
                        iright = min(level-4,nr-1)

                        mpi_ptn=0
                        do ii = proc_irange(1,myid+1,level-5),proc_irange(2,myid+1,level-5)
                        
                            jleft = max(2,level-ii-np+1)
                            jright = min(level-ii-2,nt-1)
                            do jj = jleft,jright

                                kk = level-ii-jj
                                
                                if(rdirec<0) then
                                    iir=ii
                                else
                                    iir=nr+1-ii
                                end if
                                if(tdirec<0) then
                                    iit=jj
                                else
                                    iit=nt+1-jj
                                end if
                                if(pdirec<0) then
                                    iip=kk
                                else
                                    iip=np+1-kk
                                end if
            

                                if(ischange(iir,iit,iip)==1) then
                                    

                                    !sigr=2*sqrt(a(iir,iit,iip)); sigt=2*sqrt(b(iir,iit,iip)); sigp=2*sqrt(c(iir,iit,iip))
                                    !coe=1.0/((sigr/dr)+(sigt/dt)+(sigp/dp))

                                    sigr = sqrt(a(iir,iit,iip))
                                    sigt = sqrt(b(iir,iit,iip))
                                    sigp = sqrt(c(iir,iit,iip))
                                    coe=1.0/((sigr/dr)+(sigt/dt)+(sigp/dp))

                                    ! partial derivatives

                                    pr1=(tau(iir,iit,iip)-tau(iir-1,iit,iip))/dr;
                                    pr2=(tau(iir+1,iit,iip)-tau(iir,iit,iip))/dr;                                       
                    
                                    pt1=(tau(iir,iit,iip)-tau(iir,iit-1,iip))/dt; 
                                    pt2=(tau(iir,iit+1,iip)-tau(iir,iit,iip))/dt;  
                                
                                    pp1=(tau(iir,iit,iip)-tau(iir,iit,iip-1))/dp;                    
                                    pp2=(tau(iir,iit,iip+1)-tau(iir,iit,iip))/dp;

                                    !计算 LF Hamiltonian 

                                    Htau=sqrt( a(iir,iit,iip)*((pr1+pr2)/2)**2 + b(iir,iit,iip)*((pt1+pt2)/2)**2 &
                                    & + c(iir,iit,iip)*((pp1+pp2)/2)**2 -2*f(iir,iit,iip)*(pt1+pt2)/2*(pp1+pp2)/2)
                                    
                                    ! 更新 update timetable
                                    tpT=coe*(fun(iir,iit,iip)-Htau)  &
                                    & +coe*(sigr*(pr2-pr1)/2+sigt*(pt2-pt1)/2+sigp*(pp2-pp1)/2)+tau(iir,iit,iip);

                                    tau(iir,iit,iip) = tpT

                                    
                                    mpi_ptn=mpi_ptn+1
                                    mpi_ptijk(:,mpi_ptn) = (/iir,iit,iip/)
                                    mpi_ptv(mpi_ptn) = tpT
                                    my_npt = my_npt + 1
                                end if                              
                            end do
                        end do
                        
                        !if (1>2) then
                        ! 数据通信 步骤1, 其他线程传输给主线程
                        if (myid .eq. 0) then
                            ! 负责接受数据
                            do iproc=1,nproc-1
                                call mpi_recv(int_temp,1,mpi_integer,iproc,tag,mpi_comm_world,istat,ierr)
                                if (int_temp .ne. 0) then
                                    call mpi_recv(mpi_ptijk(1,mpi_ptn+1),int_temp*3,mpi_integer,iproc,tag+1, &
                                                & mpi_comm_world,istat,ierr)
                                    call mpi_recv(mpi_ptv(mpi_ptn+1),int_temp,mpi_double_precision,iproc,tag+2, &
                                                & mpi_comm_world,istat,ierr)
                                    ! 接收完数据，赋值
                                    do ipt=mpi_ptn+1,mpi_ptn+int_temp
                                        tau(mpi_ptijk(1,ipt),mpi_ptijk(2,ipt),mpi_ptijk(3,ipt)) = mpi_ptv(ipt)
                                    end do
                                    mpi_ptn = mpi_ptn + int_temp
                                end if
                            end do
                        else
                            ! 负责发送数据
                            call mpi_send(mpi_ptn,1,mpi_integer,0,tag,mpi_comm_world,ierr)
                            !print *, mpi_ptn, mpi_ptijk(:,1:mpi_ptn), mpi_ptv(1:mpi_ptn)
                            if (mpi_ptn .ne. 0) then
                                call mpi_send(mpi_ptijk,3*mpi_ptn,mpi_integer,0,tag+1,mpi_comm_world,ierr)
                                call mpi_send(mpi_ptv,mpi_ptn,mpi_double_precision,0,tag+2,mpi_comm_world,ierr)    
                            end if                       
                        end if
                        
                        ! 数据通信 步骤2, 主线程将更新后的level数据广播到其他线程

                        call mpi_bcast(mpi_ptn,1,mpi_integer,0,mpi_comm_world,ierr)
                        call mpi_bcast(mpi_ptijk,3*mpi_ptn,mpi_integer,0,mpi_comm_world,ierr)
                        call mpi_bcast(mpi_ptv,mpi_ptn,mpi_double_precision,0,mpi_comm_world,ierr)    
                        
                        ! 数据广播完成，开始赋值
                        if (myid .ne. 0) then
                            do ipt=1,mpi_ptn
                                tau(mpi_ptijk(1,ipt),mpi_ptijk(2,ipt),mpi_ptijk(3,ipt)) = mpi_ptv(ipt)
                            end do
                        end if

                        call mpi_barrier(mpi_comm_world,ierr)
                        !end if
                    end do    

                    ! print *, my_npt, Nptn
                    
                    ! 处理边界

                    do iit=2,nt-1
                        do iip=2,np-1
                            if (tau(3,iit,iip)>0) then
                                tau(1,iit,iip) = max(2*tau(2,iit,iip)-tau(3,iit,iip),tau(3,iit,iip))
                            else
                                tau(1,iit,iip) = min(2*tau(2,iit,iip)-tau(3,iit,iip),tau(3,iit,iip))
                            end if
                            if (tau(nr-2,iit,iip)>0) then
                                tau(nr,iit,iip) = max(2*tau(nr-1,iit,iip)-tau(nr-2,iit,iip),tau(nr-2,iit,iip))
                            else
                                tau(nr,iit,iip) = min(2*tau(nr-1,iit,iip)-tau(nr-2,iit,iip),tau(nr-2,iit,iip))
                            end if
                        end do
                    end do
                    do iir=2,nr-1
                        do iip=2,np-1
                            if (tau(iir,3,iip)>0) then
                                tau(iir,1,iip) = max(2*tau(iir,2,iip)-tau(iir,3,iip),tau(iir,3,iip))
                            else
                                tau(iir,1,iip) = min(2*tau(iir,2,iip)-tau(iir,3,iip),tau(iir,3,iip))
                            end if
                            if (tau(iir,nt-2,iip)>0) then
                                tau(iir,nt,iip) = max(2*tau(iir,nt-1,iip)-tau(iir,nt-2,iip),tau(iir,nt-2,iip))
                            else
                                tau(iir,nt,iip) = min(2*tau(iir,nt-1,iip)-tau(iir,nt-2,iip),tau(iir,nt-2,iip))
                            end if
                        end do
                    end do
                    do iir=2,nr-1
                        do iit=2,nt-1
                            if (tau(iir,iit,3)>0) then
                                tau(iir,iit,1) = max(2*tau(iir,iit,2)-tau(iir,iit,3),tau(iir,iit,3))
                            else 
                                tau(iir,iit,1) = min(2*tau(iir,iit,2)-tau(iir,iit,3),tau(iir,iit,3))
                            end if
                            if (tau(iir,iit,np-2)>0) then
                                tau(iir,iit,np) = max(2*tau(iir,iit,np-1)-tau(iir,iit,np-2),tau(iir,iit,np-2))
                            else
                                tau(iir,iit,np) = min(2*tau(iir,iit,np-1)-tau(iir,iit,np-2),tau(iir,iit,np-2))
                            end if
                        end do
                    end do   


                end do
            end do
        end do

        

        ! 统计误差，判断迭代终止条件
        do iir=2,nr-1
            do iit=2,nt-1
                do iip=2,np-1
                    L1_dif=L1_dif+abs(tau(iir,iit,iip)-tau_old(iir,iit,iip))
                    Linf_dif=max(Linf_dif,abs(tau(iir,iit,iip)-tau_old(iir,iit,iip)))
                end do               
            end do
        end do
        
        do iir=3,nr-2
            do iit=3,nt-2
                do iip=3,np-2
                    L1_err=L1_err+abs(u(iir,iit,iip)-tau(iir,iit,iip))
                    Linf_err=max(Linf_err,abs(tau(iir,iit,iip)-u(iir,iit,iip)))   
                end do            
            end do
        end do
        L1_err=L1_err/((nr-4)*(nt-4)*(np-4))
        L1_dif=L1_dif/((nr-2)*(nt-2)*(np-2))

        if (myid .eq. 0) then 
            ! ################  iteration information  #################
            if (abs(L1_dif)<tol) then
                write(*,*) 'id',myid,'iter ',iter,', T is steadt'
                isexit = .true.
            else
                write(*,*) 'id',myid,'iter ',iter,', T is changing, L1 dif = ', L1_dif,'L inf dif = ', Linf_dif
            end if

            if (iter==MaxIter) then    
                write(*,*) 'id','iter ',iter,', max iteration steps'
            end if

            ! ################  solver accuracy  #################
            !write(*,'(a,f15.7,a,f15.7)') 'L_1(T-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err
        end if

        call mpi_bcast(isexit,1,mpi_logical,0,mpi_comm_world,ierr)

        if (isexit) then
            exit
        end if
    end do

    T=tau
    
end subroutine

!  3-D first order LF scheme 
subroutine FSM_O1_sphe_3d(rr,tt,pp,nr,nt,np,spha,sphb,sphc,sphf,T,fun,r0,t0,p0,u)
    integer nr,nt,np
    double precision :: dr,dt,dp,rr(nr),tt(nt),pp(np),a(nr,nt,np),b(nr,nt,np),c(nr,nt,np),f(nr,nt,np)
    double precision :: spha(nr,nt,np),sphb(nr,nt,np),sphc(nr,nt,np),sphf(nr,nt,np)
    double precision :: fun(nr,nt,np),T(nr,nt,np),ischange(nr,nt,np),u(nr,nt,np),xi(nr,nt,np),eta(nr,nt,np)
    double precision :: r1,r2,r3,a0,b0,c0,d0,f0,T0v(nr,nt,np),T0r(nr,nt,np),T0t(nr,nt,np),T0p(nr,nt,np)
    double precision :: tau(nr,nt,np),tau_old(nr,nt,np),r0,t0,p0,xi0,eta0
    integer :: idr0,idt0,idp0
    double precision :: sigr,sigt,sigp,coe,pr1,pr2,wr1,wr2,pt1,pt2,wt1,wt2,pp1,pp2,wp1,wp2,tpT
    double precision :: L1_dif,Linf_dif,L1_err,Linf_err
    double precision :: t1,p1
    double precision,parameter :: tol = (10.0)**(-5),eps=10.0**(-12)
    integer,parameter :: MaxIter=1000
    integer iter,rdirec,tdirec,pdirec,iir,iit,iip
    double precision :: x0,y0,z0,x,y,z,xst,yst,zst,e11,e12,e13,e21,e22,e23,e31,e32,e33
    double precision :: xstr,xstt,xstp,ystr,ystt,ystp,zstr,zstt,zstp

    ! ------------------------ 构造网格 build mesh ------------------------ 
    dr=rr(2)-rr(1); dt=tt(2)-tt(1); dp=pp(2)-pp(1)

    ! ------------------------ 构造矩阵 build eikonal matrix -------------------------  
    !    a -d -e
    !   -d  b -f
    !   -e -f  c
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                a(iir,iit,iip) = spha(iir,iit,iip)
                b(iir,iit,iip) = sphb(iir,iit,iip)/(rr(iir)**2)
                c(iir,iit,iip) = sphc(iir,iit,iip)/(rr(iir)**2*cos(tt(iit))**2)
                f(iir,iit,iip) = sphf(iir,iit,iip)/(rr(iir)**2*cos(tt(iit)))
            end do
        end do
    end do
    
    ! ------------------------ 构造 T0 ------------------------ 

    ! 震源处参数离散化 source parameterization
    !idr0=floor((r0-rr(1))/dr+1); idt0=floor((t0-tt(1))/dt+1); idp0=floor((p0-pp(1))/dp+1); 
    !r1 = min(1.0,(r0-rr(idr0))/dr); r2 = min(1.0,(t0-tt(idt0))/dt); r3 = min(1.0,(p0-pp(idp0))/dp); 

    idr0 = int((r0-rr(1))/dr+1); idt0=int((t0-tt(1))/dt+1); idp0=int((p0-pp(1))/dp+1); 

    ! initial tau
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                !if ( iir==idr0 .and. iit==idt0 .and.iip==idp0 ) then
                if ( ((rr(iir)-r0)/dr)**2 + ((tt(iit)-t0)/dt)**2 + ((pp(iip)-p0)/dp)**2 < 2**2) then
                    tau(iir,iit,iip) = u(iir,iit,iip)    ! source condition
                    ischange(iir,iit,iip)=0
                    if (iir==1 .or. iir==nr .or. iit==1 .or. iit==nt .or. iip==1 .or. iip==np) then
                        write(*,*) 'source on the boundary, mesh error'
                        print *, rr(iir),r0,tt(iit),t0,pp(iip),p0
                        pause
                    end if
                else
                    tau(iir,iit,iip) = 300
                    ischange(iir,iit,iip)=1
                end if
            end do
        end do
    end do


    ! 正式迭代，更新 iteration start, update tau
    do iter =1,MaxIter
        tau_old = tau
        L1_dif=0; Linf_dif=0;L1_err=0;Linf_err=0;
        do rdirec = -1,1,2
            do tdirec = -1,1,2
                do pdirec = -1,1,2

                    !x: nr-1 <-> 2, y: nt-1 <-> 2, z: np-1 <-> 2

                    do iir=nint(0.5+nr/2.0+(nr/2.0-1.5)*rdirec),nint(0.5+nr/2.0+(-nr/2.0+1.5)*rdirec),-rdirec
                        do iit=nint(0.5+nt/2.0+(nt/2.0-1.5)*tdirec),nint(0.5+nt/2.0+(-nt/2.0+1.5)*tdirec),-tdirec
                            do iip=nint(0.5+np/2.0+(np/2.0-1.5)*pdirec),nint(0.5+np/2.0+(-np/2.0+1.5)*pdirec),-pdirec
                                
                                if(ischange(iir,iit,iip)==1) then
                                    sigr=sqrt(a(iir,iit,iip)); sigt=sqrt(b(iir,iit,iip)); sigp=sqrt(c(iir,iit,iip))
                                    coe=1.0/((sigr/dr)+(sigt/dt)+(sigp/dp))

                                    ! forward and backward partial derivatives

                                    
                                        pr1=(tau(iir,iit,iip)-tau(iir-1,iit,iip))/dr;
                                        pr2=(tau(iir+1,iit,iip)-tau(iir,iit,iip))/dr;                                       
                       
                                        pt1=(tau(iir,iit,iip)-tau(iir,iit-1,iip))/dt; 
                                        pt2=(tau(iir,iit+1,iip)-tau(iir,iit,iip))/dt;  
                                   
                                        pp1=(tau(iir,iit,iip)-tau(iir,iit,iip-1))/dp;                    
                                        pp2=(tau(iir,iit,iip+1)-tau(iir,iit,iip))/dp; 


                                    ! calculate  LF Hamiltonian 

                                    Htau=sqrt( a(iir,iit,iip)*((pr1+pr2)/2)**2 + b(iir,iit,iip)*((pt1+pt2)/2)**2 &
                                      & + c(iir,iit,iip)*((pp1+pp2)/2)**2 -2*f(iir,iit,iip)*(pt1+pt2)/2*(pp1+pp2)/2 )
                                    
                                    ! 更新 update timetable
                                    tpT=coe*(fun(iir,iit,iip)-Htau)  &
                                     & +coe*(sigr*(pr2-pr1)/2+sigt*(pt2-pt1)/2+sigp*(pp2-pp1)/2)+tau(iir,iit,iip);

                                    !write(*,*) fun(iir,iit,iip),Htau

                                     if (tpT < tau(iir,iit,iip)) then
                                        tau(iir,iit,iip) = tpT
                                     end if

                                end if

                                
                            end do
                        end do
                    end do    

                    ! boundary

                    do iit=1,nt
                        do iip=1,np
                            tau(1,iit,iip) = max(2*tau(2,iit,iip)-tau(3,iit,iip),tau(3,iit,iip))
                            tau(nr,iit,iip) = max(2*tau(nr-1,iit,iip)-tau(nr-2,iit,iip),tau(nr-2,iit,iip))
                        end do
                    end do
                    do iir=1,nr
                        do iip=1,np
                            tau(iir,1,iip) = max(2*tau(iir,2,iip)-tau(iir,3,iip),tau(iir,3,iip))
                            tau(iir,nt,iip) = max(2*tau(iir,nt-1,iip)-tau(iir,nt-2,iip),tau(iir,nt-2,iip))
                        end do
                    end do
                    do iir=1,nr
                        do iit=1,nt
                            tau(iir,iit,1) = max(2*tau(iir,iit,2)-tau(iir,iit,3),tau(iir,iit,3))
                            tau(iir,iit,np) = max(2*tau(iir,iit,np-1)-tau(iir,iit,np-2),tau(iir,iit,np-2))
                        end do
                    end do   


                end do
            end do
        end do

        

        ! 统计误差，判断迭代终止条件
        do iir=2,nr-1
            do iit=2,nt-1
                do iip=2,np-1
                    L1_dif=L1_dif+abs(tau(iir,iit,iip)-tau_old(iir,iit,iip))
                    Linf_dif=max(Linf_dif,abs(tau(iir,iit,iip)-tau_old(iir,iit,iip)))
                end do               
            end do
        end do
        
        do iir=3,nr-2
            do iit=3,nt-2
                do iip=3,np-2
                    L1_err=L1_err+abs(u(iir,iit,iip)-tau(iir,iit,iip))
                    Linf_err=max(Linf_err,abs(tau(iir,iit,iip)-u(iir,iit,iip)))   
                end do            
            end do
        end do
        L1_err=L1_err/((nr-4)*(nt-4)*(np-4))
        L1_dif=L1_dif/((nr-2)*(nt-2)*(np-2))

        if (myid .eq. 0) then 
            ! ################  iteration information  #################
            if (abs(L1_dif)<tol) then
                write(*,*) 'iter ',iter,', T is steadt'
                exit
            else
                write(*,*) 'iter ',iter,', T is changing, L1 dif = ', L1_dif,'L inf dif = ', Linf_dif
            end if

            if (iter==MaxIter) then    
                write(*,*) 'iter ',iter,', max iteration steps'
            end if

            ! ################  solver accuracy  #################
            !write(*,'(a,f15.7,a,f15.7)') 'L_1(T-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err
        end if

    end do

    T=tau
    
end subroutine

!  3-D interpolation non-uniform mesh
subroutine Linear_Interp_3D_nonuni(xx,yy,zz,val,nx,ny,nz,x0,y0,z0,v0)
    integer :: nx,ny,nz,idx0,idy0,idz0
    double precision :: xx(nx),yy(ny),zz(nz),val(nx,ny,nz),x0,y0,z0,v0,dx,dy,dz
    double precision :: r1,r2,r3
    integer :: iix,iiy,iiz

    idx0=-1; idy0=-1; idz0=-1
    do iix=1,nx-1
        if (xx(iix)<= x0 .and. xx(iix+1)>x0 ) then
            idx0=iix
            exit
        end if
    end do
    do iiy=1,ny-1
        if (yy(iiy)<= y0 .and. yy(iiy+1)>y0 ) then
            idy0=iiy
            exit
        end if
    end do
    do iiy=1,ny-1
        if (zz(iiz)<= z0 .and. zz(iiz+1)>z0 ) then
            idz0=iiz
            exit
        end if
    end do
    
    if (idx0<=0 .or. idy0<0 .or. idz0<0) then
        write(*,*) 'point out of the mesh'
        pause
        return
    end if

    dx = xx(idx0+1)-xx(idx0)
    dy = yy(idy0+1)-yy(idy0)
    dz = zz(idz0+1)-zz(idz0)
    
    r1 = min(1.0, (x0-xx(idx0))/dx )
    r2 = min(1.0, (y0-yy(idy0))/dy )
    r3 = min(1.0, (z0-zz(idz0))/dz )

    v0 = (1-r1)*(1-r2)*(1-r3)*val(idx0,idy0,idz0) + (1-r1)*(1-r2)*r3*val(idx0,idy0,idz0+1) &
     & + (1-r1)*r2*(1-r3)*val(idx0,idy0+1,idz0) + (1-r1)*r2*r3*val(idx0,idy0+1,idz0+1) &
     & + r1*(1-r2)*(1-r3)*val(idx0+1,idy0,idz0) + r1*(1-r2)*r3*val(idx0+1,idy0,idz0+1) &
     & + r1*r2*(1-r3)*val(idx0+1,idy0+1,idz0) + r1*r2*r3*val(idx0+1,idy0+1,idz0+1)
end subroutine

!  3-D interpolation uniform mesh
subroutine Linear_Interp_3D(xx,yy,zz,val,nx,ny,nz,x0,y0,z0,v0)
    integer :: nx,ny,nz,idx0,idy0,idz0
    double precision :: xx(nx),yy(ny),zz(nz),val(nx,ny,nz),x0,y0,z0,v0,dx,dy,dz
    double precision :: r1,r2,r3
    integer :: iix,iiy,iiz
    
    dx = xx(2)-xx(1)
    dy = yy(2)-yy(1)
    dz = zz(2)-zz(1)

    idx0 = floor((x0-xx(1))/dx)+1
    idy0 = floor((y0-yy(1))/dy)+1
    idz0 = floor((z0-zz(1))/dz)+1

    if (idx0<1 .or. idx0>nx-1 .or. idy0<1 .or. idy0>ny-1 .or. idz0<1 .or. idz0>nz-1) then
        write(*,*) 'point out of the mesh'
        pause
        return
    end if

    r1 = min(1.0, (x0-xx(idx0))/dx )
    r2 = min(1.0, (y0-yy(idy0))/dy )
    r3 = min(1.0, (z0-zz(idz0))/dz )

    v0 = (1-r1)*(1-r2)*(1-r3)*val(idx0,idy0,idz0) + (1-r1)*(1-r2)*r3*val(idx0,idy0,idz0+1) &
        & + (1-r1)*r2*(1-r3)*val(idx0,idy0+1,idz0) + (1-r1)*r2*r3*val(idx0,idy0+1,idz0+1) &
        & + r1*(1-r2)*(1-r3)*val(idx0+1,idy0,idz0) + r1*(1-r2)*r3*val(idx0+1,idy0,idz0+1) &
        & + r1*r2*(1-r3)*val(idx0+1,idy0+1,idz0) + r1*r2*r3*val(idx0+1,idy0+1,idz0+1)

end subroutine

subroutine FSM_O1_Adj_sphe_3d(rr,tt,pp,nr,nt,np,Table,TableADJ,zeta,xi,eta,Nrec,rrec,trec,prec,sourceADJ)

    integer :: nr,nt,np,Nrec
    double precision :: dr,dt,dp,rr(nr),tt(nt),pp(np)
    double precision :: Table(nr,nt,np),TableADJ(nr,nt,np),delta(nr,nt,np)
    double precision :: zeta(nr,nt,np),xi(nr,nt,np),eta(nr,nt,np)
    double precision :: rrec(Nrec),trec(Nrec),prec(Nrec),sourceADJ(Nrec)
    double precision :: a1,a1m(nr,nt,np),a1p(nr,nt,np),a2,a2m(nr,nt,np),a2p(nr,nt,np)
    double precision :: b1,b1m(nr,nt,np),b1p(nr,nt,np),b2,b2m(nr,nt,np),b2p(nr,nt,np)
    double precision :: c1,c1m(nr,nt,np),c1p(nr,nt,np),c2,c2m(nr,nt,np),c2p(nr,nt,np)
    double precision :: coe

    integer,parameter :: MaxIter=100
    double precision,parameter :: tol=10.0**(-6),eps=10.0**(-6)

    integer :: iir,iit,iip,idi,idj,idk,rdirec,tdirec,pdirec
    double precision :: r1,r2,r3,Linf_dif,tpTabldADJ,Hamilton
    double precision :: tmpr1,tmpr2,tmpt1,tmpt2

    ! #######  网格间距  mesh size  ########
    dr = rr(2)-rr(1) 
    dt = tt(2)-tt(1) 
    dp = pp(2)-pp(1)
    
    ! #################  构造源项 build the sourece term, the delta function  #################
    do iir = 1,nr
        do iit = 1,nt
            do iip= 1,np
                delta(iir,iit,iip) = 0
                TableADJ(iir,iit,iip) = 0
            end do
        end do
    end do

    ! ------- 遍历每个台站，给予delta函数贡献 loop each station to contribute the delta function -------
    do ir=1,Nrec  
        idi=floor((rrec(ir)-rr(1))/dr+1); 
        idj=floor((trec(ir)-tt(1))/dt+1);
        idk=floor((prec(ir)-pp(1))/dp+1);

        r1 = min(1.0,(rrec(ir)-rr(idi))/dr); 
        r2 = min(1.0,(trec(ir)-tt(idj))/dt); 
        r3 = min(1.0,(prec(ir)-pp(idk))/dp); 

        delta(idi,idj,idk) = delta(idi,idj,idk) + sourceADJ(ir)*(1-r1)*(1-r2)*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj,idk) = delta(idi+1,idj,idk) + sourceADJ(ir)*r1*(1-r2)*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi,idj+1,idk) = delta(idi,idj+1,idk) + sourceADJ(ir)*(1-r1)*r2*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj+1,idk) = delta(idi+1,idj+1,idk) + sourceADJ(ir)*r1*r2*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi,idj,idk+1) = delta(idi,idj,idk+1) + sourceADJ(ir)*(1-r1)*(1-r2)*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj,idk+1) = delta(idi+1,idj,idk+1) + sourceADJ(ir)*r1*(1-r2)*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi,idj+1,idk+1) = delta(idi,idj+1,idk+1) + sourceADJ(ir)*(1-r1)*r2*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj+1,idk+1) = delta(idi+1,idj+1,idk+1) + sourceADJ(ir)*r1*r2*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        !print *, delta(idi,idj,idk),delta(idi+1,idj+1,idk+1)
    end do

    ! #################  构造方程系数 calculate coefficients of equations  #################
    do iir = 2,nr-1
        do iit = 2,nt-1
            do iip= 2,np-1
                tmpr1 = (rr(iir-1)+rr(iir))/2; tmpt1 = (tt(iit-1)+tt(iit))/2;
                tmpr2 = (rr(iir)+rr(iir+1))/2; tmpt2 = (tt(iit)+tt(iit-1))/2;

                a1 = -(1+zeta(iir-1,iit,iip)+zeta(iir,iit,iip))*(Table(iir,iit,iip)-Table(iir-1,iit,iip))/dr
                a1m(iir,iit,iip) = (a1-abs(a1))/2; a1p(iir,iit,iip) = (a1+abs(a1))/2;
                a2 = -(1+zeta(iir,iit,iip)+zeta(iir+1,iit,iip))*(Table(iir+1,iit,iip)-Table(iir,iit,iip))/dr
                a2m(iir,iit,iip) = (a2-abs(a2))/2; a2p(iir,iit,iip) = (a2+abs(a2))/2;

                !b1 = -(1-xi(iir,iit-1,iip)-xi(iir,iit,iip))/(tmpr1**2)*(Table(iir,iit,iip)-Table(iir,iit-1,iip))/dt &
                !   & -(eta(iir,iit-1,iip)+eta(iir,iit,iip))/(cos(tmpt1)*tmpr1**2)/(4*dp) &
                b1 = -(1-xi(iir,iit-1,iip)-xi(iir,iit,iip))/(rr(iir)**2)*(Table(iir,iit,iip)-Table(iir,iit-1,iip))/dt &
                   & -(eta(iir,iit-1,iip)+eta(iir,iit,iip))/(cos(tmpt1)*rr(iir)**2)/(4*dp) &
                   & *((Table(iir,iit-1,iip+1)-Table(iir,iit-1,iip-1))+(Table(iir,iit,iip+1)-Table(iir,iit,iip-1)))
                b1m(iir,iit,iip) = (b1-abs(b1))/2; b1p(iir,iit,iip) = (b1+abs(b1))/2;

                !b2 = -(1-xi(iir,iit,iip)-xi(iir,iit+1,iip))/(tmpr2**2)*(Table(iir,iit+1,iip)-Table(iir,iit,iip))/dt &
                !   & -(eta(iir,iit,iip)+eta(iir,iit+1,iip))/(cos(tmpt2)*tmpr2**2)/(4*dp) &
                b2 = -(1-xi(iir,iit,iip)-xi(iir,iit+1,iip))/(rr(iir)**2)*(Table(iir,iit+1,iip)-Table(iir,iit,iip))/dt &
                   & -(eta(iir,iit,iip)+eta(iir,iit+1,iip))/(cos(tmpt2)*rr(iir)**2)/(4*dp) &
                   & *((Table(iir,iit,iip+1)-Table(iir,iit,iip-1))+(Table(iir,iit+1,iip+1)-Table(iir,iit+1,iip-1)))
                b2m(iir,iit,iip) = (b2-abs(b2))/2; b2p(iir,iit,iip) = (b2+abs(b2))/2;

                !c1 = -(eta(iir,iit,iip-1)+eta(iir,iit,iip))/(cos(tmpt1)*tmpr1**2)/(4*dt) &
                c1 = -(eta(iir,iit,iip-1)+eta(iir,iit,iip))/(cos(tt(iit))*rr(iir)**2)/(4*dt) &
                   & *((Table(iir,iit+1,iip-1)-Table(iir,iit-1,iip-1))+(Table(iir,iit+1,iip)-Table(iir,iit-1,iip))) &
                !   & -(1+xi(iir,iit,iip-1)+xi(iir,iit,iip))/(cos(tmpt1)**2*tmpr1**2)*(Table(iir,iit,iip)-Table(iir,iit,iip-1))/dp
                   & -(1+xi(iir,iit,iip-1)+xi(iir,iit,iip))/ &
                   & (cos(tt(iit))**2*rr(iir)**2)*(Table(iir,iit,iip)-Table(iir,iit,iip-1))/dp
                c1m(iir,iit,iip) = (c1-abs(c1))/2; c1p(iir,iit,iip) = (c1+abs(c1))/2;

                !c2 = -(eta(iir,iit,iip)+eta(iir,iit,iip+1))/(cos(tmpt2)*tmpr2**2)/(4*dt) &
                c2 = -(eta(iir,iit,iip)+eta(iir,iit,iip+1))/(cos(tt(iit))*rr(iir)**2)/(4*dt) &
                   & *((Table(iir,iit+1,iip)-Table(iir,iit-1,iip))+(Table(iir,iit+1,iip+1)-Table(iir,iit-1,iip+1))) &
                !   & -(1+xi(iir,iit,iip)+xi(iir,iit,iip+1))/(cos(tmpt2)**2*tmpr2**2)*(Table(iir,iit,iip+1)-Table(iir,iit,iip))/dp
                   & -(1+xi(iir,iit,iip)+xi(iir,iit,iip+1)) &
                   & /(cos(tt(iit))**2*rr(iir)**2)*(Table(iir,iit,iip+1)-Table(iir,iit,iip))/dp
                c2m(iir,iit,iip) = (c2-abs(c2))/2; c2p(iir,iit,iip) = (c2+abs(c2))/2;
                
            end do
        end do
    end do

    
    ! ################# 固定伴随场边界条件 fix the boundary of the adjoint field ######################
    ! -------- R,Theta boundary----- 
    do iir=1,nr
        do iit=1,nt
            TableADJ(iir,iit,1) = 0
            TableADJ(iir,iit,np) = 0
        end do
    end do
    ! -------- R,Phi boundary----- 
    do iir=1,nr
        do iit=1,nt
            TableADJ(iir,1,iip) = 0
            TableADJ(iir,nt,iip) = 0
        end do
    end do
    ! -------- Theta,Phi boundary----- 
    do iir=1,nr
        do iit=1,nt
            TableADJ(1,iit,iip) = 0
            TableADJ(nr,iit,iip) = 0
        end do
    end do




    ! ################ 使用FSM计算伴随场 calculate adjoint field by FSM ######################

    do iter =1,MaxIter
        Linf_dif=0;
        do rdirec = -1,1,2
            do tdirec = -1,1,2
                do pdirec = -1,1,2

                    !x: nr-1 <-> 2, y: nt-1 <-> 2, z: np-1 <-> 2
                    do iir=nint(0.5+nr/2.0+(nr/2.0-1.5)*rdirec),nint(0.5+nr/2.0+(-nr/2.0+1.5)*rdirec),-rdirec
                        do iit=nint(0.5+nt/2.0+(nt/2.0-1.5)*tdirec),nint(0.5+nt/2.0+(-nt/2.0+1.5)*tdirec),-tdirec
                            do iip=nint(0.5+np/2.0+(np/2.0-1.5)*pdirec),nint(0.5+np/2.0+(-np/2.0+1.5)*pdirec),-pdirec
                                
                                
                                coe = (a2p(iir,iit,iip)-a1m(iir,iit,iip))/dr &
                                   & +(b2p(iir,iit,iip)-b1m(iir,iit,iip))/dt &
                                   & +(c2p(iir,iit,iip)-c1m(iir,iit,iip))/dp
                                 
                                if (abs(coe)<eps) then
                                    TableADJ(iir,iit,iip) = 0
                                else
                                    Hamilton = &
                                    &   (TableADJ(iir-1,iit,iip)*a1p(iir,iit,iip)-TableADJ(iir+1,iit,iip)*a2m(iir,iit,iip))/dr &
                                    & + (TableADJ(iir,iit-1,iip)*b1p(iir,iit,iip)-TableADJ(iir,iit+1,iip)*b2m(iir,iit,iip))/dt &
                                    & + (TableADJ(iir,iit,iip-1)*c1p(iir,iit,iip)-TableADJ(iir,iit,iip+1)*c2m(iir,iit,iip))/dp

                                    tpTableADJ = (delta(iir,iit,iip)+Hamilton)/coe
                                    

                                    !if (abs(tpTableADJ-TableADJ(iir,iit,iip))>Linf_dif) then
                                    !    print *, iir,iit,iip,tpTableADJ,TableADJ(iir,iit,iip)
                                    !end if
                                    Linf_dif = max(Linf_dif, abs(tpTableADJ-TableADJ(iir,iit,iip)))
                                    
                                    TableADJ(iir,iit,iip) = tpTableADJ
                                end if
                                
                            end do
                        end do
                    end do 

                end do
            end do
        end do

        

        if (abs(Linf_dif)<tol) then
            write(*,*) 'iter ',iter,', TableADJ is steady'
            exit
        else
            write(*,*) 'iter ',iter,', TableADJ is changing, continue ... '
        end if

        if (iter==MaxIter) then    
        !    write(*,*) 'iter ',iter,', max iteration steps'
        end if

    end do
    
end subroutine

subroutine FSM_O1_Adj_PS_sphe_3d(rr,tt,pp,nr,nt,np,T0para,tau,TableADJ,zeta,xi,eta,Nrec,rrec,trec,prec,sourceADJ)

    integer :: nr,nt,np,Nrec
    double precision :: dr,dt,dp,rr(nr),tt(nt),pp(np)
    double precision :: tau(nr,nt,np),Table(nr,nt,np),TableADJ(nr,nt,np),delta(nr,nt,np)
    double precision :: zeta(nr,nt,np),xi(nr,nt,np),eta(nr,nt,np)
    double precision :: rrec(Nrec),trec(Nrec),prec(Nrec),sourceADJ(Nrec)
    double precision :: a1,a1m(nr,nt,np),a1p(nr,nt,np),a2,a2m(nr,nt,np),a2p(nr,nt,np)
    double precision :: b1,b1m(nr,nt,np),b1p(nr,nt,np),b2,b2m(nr,nt,np),b2p(nr,nt,np)
    double precision :: c1,c1m(nr,nt,np),c1p(nr,nt,np),c2,c2m(nr,nt,np),c2p(nr,nt,np)
    double precision :: coe

    integer,parameter :: MaxIter=100
    double precision,parameter :: tol=10.0**(-6),eps=10.0**(-6)

    integer :: iir,iit,iip,idi,idj,idk,rdirec,tdirec,pdirec
    double precision :: r1,r2,r3,Linf_dif,tpTabldADJ,Hamilton
    double precision :: tmpr1,tmpr2,tmpt1,tmpt2

    double precision :: r0,t0,p0,t1,p1
    double precision :: T0para(8),a0,b0,c0,f0,fun0
    double precision :: T0v_hr(nr-1,nt,np),T0r_hr(nr-1,nt,np),T0v(nr,nt,np)
    double precision :: T0v_ht(nr,nt-1,np),T0t_ht(nr,nt-1,np),T0p_ht(nr,nt-1,np)
    double precision :: T0v_hp(nr,nt,np-1),T0t_hp(nr,nt,np-1),T0p_hp(nr,nt,np-1)

    ! #######  网格间距  mesh size  ########
    dr = rr(2)-rr(1) 
    dt = tt(2)-tt(1) 
    dp = pp(2)-pp(1)
    
    ! #################  构造源项 build the sourece term, the delta function  #################
    do iir = 1,nr
        do iit = 1,nt
            do iip= 1,np
                delta(iir,iit,iip) = 0
            end do
        end do
    end do

    ! ------- 遍历每个台站，给予delta函数贡献 loop each station to contribute the delta function -------
    do ir=1,Nrec  
        idi=floor((rrec(ir)-rr(1))/dr+1); 
        idj=floor((trec(ir)-tt(1))/dt+1);
        idk=floor((prec(ir)-pp(1))/dp+1);

        r1 = min(1.0,(rrec(ir)-rr(idi))/dr); 
        r2 = min(1.0,(trec(ir)-tt(idj))/dt); 
        r3 = min(1.0,(prec(ir)-pp(idk))/dp); 

        delta(idi,idj,idk) = delta(idi,idj,idk) + sourceADJ(ir)*(1-r1)*(1-r2)*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj,idk) = delta(idi+1,idj,idk) + sourceADJ(ir)*r1*(1-r2)*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi,idj+1,idk) = delta(idi,idj+1,idk) + sourceADJ(ir)*(1-r1)*r2*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj+1,idk) = delta(idi+1,idj+1,idk) + sourceADJ(ir)*r1*r2*(1-r3)/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi,idj,idk+1) = delta(idi,idj,idk+1) + sourceADJ(ir)*(1-r1)*(1-r2)*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj,idk+1) = delta(idi+1,idj,idk+1) + sourceADJ(ir)*r1*(1-r2)*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi,idj+1,idk+1) = delta(idi,idj+1,idk+1) + sourceADJ(ir)*(1-r1)*r2*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        delta(idi+1,idj+1,idk+1) = delta(idi+1,idj+1,idk+1) + sourceADJ(ir)*r1*r2*r3/(dr*dp*dt*rrec(ir)**2*cos(trec(ir)))
        !print *, delta(idi,idj,idk),delta(idi+1,idj+1,idk+1)
    end do

    ! 构造T0_hr  iir,iit,iip -> r1+(iir+0.5)*dr, t1+iit*dt, p1+iip*dp  
    a0 = T0para(1); b0 = T0para(2); c0 = T0para(3); f0 = T0para(4); fun0 = T0para(5);
    r0 = T0para(6); t0 = T0para(7); p0 = T0para(8);

    do iir=1,nr-1
        do iit=1,nt
            do iip=1,np
                r1 = rr(iir)+0.5*dr; t1 = tt(iit); p1 = pp(iip)
                T0v_hr(iir,iit,iip) = fun0*sqrt( 1.0/a0*(r1-r0)**2 + c0/(c0*b0-f0**2)*(t1-t0)**2 &
                                          & + b0/(c0*b0-f0**2)*(p1-p0)**2 + 2*f0/(c0*b0-f0**2)*(t1-t0)*(p1-p0) ) 
                if ( T0v_hr(iir,iit,iip)==0 ) then
                    T0r_hr(iir,iit,iip) = 0
                else
                    T0r_hr(iir,iit,iip) = fun0**2*(1.0/a0*(r1-r0))/T0v_hr(iir,iit,iip)  
                end if
            end do
        end do
    end do
    ! 构造T0_ht  iir,iit,iip -> r1+iir*dr, t1+(iit+0.5)*dt, p1+iip*dp  
    do iir=1,nr
        do iit=1,nt-1
            do iip=1,np
                r1 = rr(iir); t1 = tt(iit)+0.5*dt; p1 = pp(iip)
                T0v_ht(iir,iit,iip) = fun0*sqrt( 1.0/a0*(r1-r0)**2 + c0/(c0*b0-f0**2)*(t1-t0)**2 &
                                          & + b0/(c0*b0-f0**2)*(p1-p0)**2 + 2*f0/(c0*b0-f0**2)*(t1-t0)*(p1-p0) ) 
                if ( T0v_ht(iir,iit,iip)==0 ) then
                    T0t_ht(iir,iit,iip) = 0
                    T0p_ht(iir,iit,iip) = 0
                else
                    T0t_ht(iir,iit,iip) = fun0**2*(c0/(c0*b0-f0**2)*(t1-t0)+f0/(c0*b0-f0**2)*(p1-p0))/T0v_ht(iir,iit,iip)
                    T0p_ht(iir,iit,iip) = fun0**2*(b0/(c0*b0-f0**2)*(p1-p0)+f0/(c0*b0-f0**2)*(t1-t0))/T0v_ht(iir,iit,iip)    
                end if
            end do
        end do
    end do
    ! 构造T0_hp  iir,iit,iip -> r1+iir*dr, t1+iit*dt, p1+(iip+0.5)*dp  
    do iir=1,nr
        do iit=1,nt
            do iip=1,np-1
                r1 = rr(iir); t1 = tt(iit); p1 = pp(iip)+0.5*dp
                T0v_hp(iir,iit,iip) = fun0*sqrt( 1.0/a0*(r1-r0)**2 + c0/(c0*b0-f0**2)*(t1-t0)**2 &
                                          & + b0/(c0*b0-f0**2)*(p1-p0)**2 + 2*f0/(c0*b0-f0**2)*(t1-t0)*(p1-p0) ) 
                if ( T0v_hp(iir,iit,iip)==0 ) then
                    T0t_hp(iir,iit,iip) = 0
                    T0p_hp(iir,iit,iip) = 0
                else
                    T0t_hp(iir,iit,iip) = fun0**2*(c0/(c0*b0-f0**2)*(t1-t0)+f0/(c0*b0-f0**2)*(p1-p0))/T0v_hp(iir,iit,iip)
                    T0p_hp(iir,iit,iip) = fun0**2*(b0/(c0*b0-f0**2)*(p1-p0)+f0/(c0*b0-f0**2)*(t1-t0))/T0v_hp(iir,iit,iip)    
                end if
            end do
        end do
    end do




    ! #################  构造方程系数 calculate coefficients of equations  #################
    do iir = 2,nr-1
        do iit = 2,nt-1
            do iip= 2,np-1
                tmpt1 = (tt(iit-1)+tt(iit))/2;
                tmpt2 = (tt(iit)+tt(iit-1))/2;

                a1 = -(1+zeta(iir-1,iit,iip)+zeta(iir,iit,iip)) &
                   & *(T0r_hr(iir-1,iit,iip)*(tau(iir,iit,iip)+tau(iir-1,iit,iip))/2 &
                   & + T0v_hr(iir-1,iit,iip)*(tau(iir,iit,iip)-tau(iir-1,iit,iip))/dr)
                a1m(iir,iit,iip) = (a1-abs(a1))/2; a1p(iir,iit,iip) = (a1+abs(a1))/2;
                a2 = -(1+zeta(iir,iit,iip)+zeta(iir+1,iit,iip)) &
                   & *(T0r_hr(iir,iit,iip)*(tau(iir+1,iit,iip)+tau(iir,iit,iip))/2 &
                   & + T0v_hr(iir,iit,iip)*(tau(iir+1,iit,iip)-tau(iir,iit,iip))/dr)
                a2m(iir,iit,iip) = (a2-abs(a2))/2; a2p(iir,iit,iip) = (a2+abs(a2))/2;

                b1 = -(1-xi(iir,iit-1,iip)-xi(iir,iit,iip))/(rr(iir)**2)* &
                   &  (T0t_ht(iir,iit-1,iip)*(tau(iir,iit,iip)+tau(iir,iit-1,iip))/2 &
                   & + T0v_ht(iir,iit-1,iip)*(tau(iir,iit,iip)-tau(iir,iit-1,iip))/dt) &
                   & -(eta(iir,iit-1,iip)+eta(iir,iit,iip))/(cos(tmpt1)*rr(iir)**2)* &
                   &  (T0p_ht(iir,iit-1,iip)*(tau(iir,iit,iip)+tau(iir,iit-1,iip))/2 &
                   & + T0v_ht(iir,iit-1,iip)*(tau(iir,iit-1,iip+1)-tau(iir,iit-1,iip-1) &
                   &                        + tau(iir,iit,iip+1)-tau(iir,iit,iip-1))/(4*dp))  
                b1m(iir,iit,iip) = (b1-abs(b1))/2; b1p(iir,iit,iip) = (b1+abs(b1))/2;
                b2 = -(1-xi(iir,iit,iip)-xi(iir,iit+1,iip))/(rr(iir)**2)* &
                   &  (T0t_ht(iir,iit,iip)*(tau(iir,iit+1,iip)+tau(iir,iit,iip))/2 &
                   & + T0v_ht(iir,iit,iip)*(tau(iir,iit+1,iip)-tau(iir,iit,iip))/dt) &
                   & -(eta(iir,iit,iip)+eta(iir,iit+1,iip))/(cos(tmpt2)*rr(iir)**2)* &
                   &  (T0p_ht(iir,iit,iip)*(tau(iir,iit+1,iip)+tau(iir,iit,iip))/2 &
                   & + T0v_ht(iir,iit,iip)*(tau(iir,iit,iip+1)-tau(iir,iit,iip-1) &
                   &                        + tau(iir,iit+1,iip+1)-tau(iir,iit+1,iip-1))/(4*dp))
                b2m(iir,iit,iip) = (b2-abs(b2))/2; b2p(iir,iit,iip) = (b2+abs(b2))/2;

                c1 = -(eta(iir,iit,iip-1)+eta(iir,iit,iip))/(cos(tt(iit))*rr(iir)**2)* &
                   &  (T0t_hp(iir,iit,iip-1)*(tau(iir,iit,iip)+tau(iir,iit,iip-1))/2 &
                   & + T0v_hp(iir,iit,iip-1)*(tau(iir,iit+1,iip-1)-tau(iir,iit-1,iip-1) &
                                           & +tau(iir,iit+1,iip)-tau(iir,iit-1,iip))/(4*dt)) &
                   & -(1+xi(iir,iit,iip-1)+xi(iir,iit,iip))/(cos(tt(iit))**2*rr(iir)**2)* &
                   &  (T0p_hp(iir,iit,iip-1)*(tau(iir,iit,iip)+tau(iir,iit,iip-1))/2 &
                   & + T0v_hp(iir,iit,iip-1)*(tau(iir,iit,iip)-tau(iir,iit,iip-1))/dp) 
                c1m(iir,iit,iip) = (c1-abs(c1))/2; c1p(iir,iit,iip) = (c1+abs(c1))/2;
                c2 = -(eta(iir,iit,iip)+eta(iir,iit,iip+1))/(cos(tt(iit))*rr(iir)**2)* &
                   &  (T0t_hp(iir,iit,iip)*(tau(iir,iit,iip+1)+tau(iir,iit,iip))/2 &
                   & + T0v_hp(iir,iit,iip)*(tau(iir,iit+1,iip)-tau(iir,iit-1,iip) &
                                         & +tau(iir,iit+1,iip+1)-tau(iir,iit-1,iip+1))/(4*dt)) &
                   & -(1+xi(iir,iit,iip)+xi(iir,iit,iip+1))/(cos(tt(iit))**2*rr(iir)**2)* &
                   &  (T0p_hp(iir,iit,iip)*(tau(iir,iit,iip+1)+tau(iir,iit,iip))/2 &
                   & + T0v_hp(iir,iit,iip)*(tau(iir,iit,iip+1)-tau(iir,iit,iip))/dp)
                c2m(iir,iit,iip) = (c2-abs(c2))/2; c2p(iir,iit,iip) = (c2+abs(c2))/2;   
                              
            end do
        end do
    end do

    
    ! ################# 固定伴随场边界条件 fix the boundary of the adjoint field ######################
    ! -------- R,Theta boundary----- 
    do iir=1,nr
        do iit=1,nt
            TableADJ(iir,iit,1) = 0
            TableADJ(iir,iit,np) = 0
        end do
    end do
    ! -------- R,Phi boundary----- 
    do iir=1,nr
        do iit=1,nt
            TableADJ(iir,1,iip) = 0
            TableADJ(iir,nt,iip) = 0
        end do
    end do
    ! -------- Theta,Phi boundary----- 
    do iir=1,nr
        do iit=1,nt
            TableADJ(1,iit,iip) = 0
            TableADJ(nr,iit,iip) = 0
        end do
    end do




    ! ################ 使用FSM计算伴随场 calculate adjoint field by FSM ######################

    do iter =1,MaxIter
        Linf_dif=0;
        do rdirec = -1,1,2
            do tdirec = -1,1,2
                do pdirec = -1,1,2

                    !x: nr-1 <-> 2, y: nt-1 <-> 2, z: np-1 <-> 2
                    do iir=nint(0.5+nr/2.0+(nr/2.0-1.5)*rdirec),nint(0.5+nr/2.0+(-nr/2.0+1.5)*rdirec),-rdirec
                        do iit=nint(0.5+nt/2.0+(nt/2.0-1.5)*tdirec),nint(0.5+nt/2.0+(-nt/2.0+1.5)*tdirec),-tdirec
                            do iip=nint(0.5+np/2.0+(np/2.0-1.5)*pdirec),nint(0.5+np/2.0+(-np/2.0+1.5)*pdirec),-pdirec
                                
                                
                                coe = (a2p(iir,iit,iip)-a1m(iir,iit,iip))/dr &
                                   & +(b2p(iir,iit,iip)-b1m(iir,iit,iip))/dt &
                                   & +(c2p(iir,iit,iip)-c1m(iir,iit,iip))/dp
                                 
                                if (abs(coe)<eps) then
                                    TableADJ(iir,iit,iip) = 0
                                else
                                    Hamilton = &
                                    &   (TableADJ(iir-1,iit,iip)*a1p(iir,iit,iip)-TableADJ(iir+1,iit,iip)*a2m(iir,iit,iip))/dr &
                                    & + (TableADJ(iir,iit-1,iip)*b1p(iir,iit,iip)-TableADJ(iir,iit+1,iip)*b2m(iir,iit,iip))/dt &
                                    & + (TableADJ(iir,iit,iip-1)*c1p(iir,iit,iip)-TableADJ(iir,iit,iip+1)*c2m(iir,iit,iip))/dp

                                    tpTableADJ = (delta(iir,iit,iip)+Hamilton)/coe
                                    

                                    !if (abs(tpTableADJ-TableADJ(iir,iit,iip))>Linf_dif) then
                                    !    print *, iir,iit,iip,tpTableADJ,TableADJ(iir,iit,iip)
                                    !end if
                                    Linf_dif = max(Linf_dif, abs(tpTableADJ-TableADJ(iir,iit,iip)))
                                    
                                    TableADJ(iir,iit,iip) = tpTableADJ
                                end if
                                
                            end do
                        end do
                    end do 

                end do
            end do
        end do

        

        if (abs(Linf_dif)<tol) then
            write(*,*) 'iter ',iter,', TableADJ is steady'
            exit
        else
        !    write(*,*) 'iter ',iter,', TableADJ is changing, continue ... '
        end if

        if (iter==MaxIter) then    
        !    write(*,*) 'iter ',iter,', max iteration steps'
        end if

    end do
    
end subroutine

subroutine Sensitivity_Kernel(rr,tt,pp,nr,nt,np,Table,TableADJ,gamma,xi,eta,fun,Ks,Kxi,Keta)
    integer :: nr,nt,np
    double precision :: rr(nr),tt(nt),pp(np),dr,dt,dp
    double precision :: Table(nr,nt,np),TableADJ(nr,nt,np)
    double precision :: gamma,xi(nr,nt,np),eta(nr,nt,np),fun(nr,nt,np)
    double precision :: Ks(nr,nt,np),Kxi(nr,nt,np),Keta(nr,nt,np)
    double precision :: Tr,Ttheta,Tphi

    integer :: iir,iit,iip

    dr = rr(2)-rr(1); dt = tt(2)-tt(1); dp = pp(2)-pp(1)
    ! inner points
    do iir=2,nr-1
        do iit=2,nt-1
            do iip=2,np-1
                Tr = (Table(iir+1,iit,iip)-Table(iir-1,iit,iip))/(2*dr)
                Ttheta = (Table(iir,iit+1,iip)-Table(iir,iit-1,iip))/(2*dt)
                Tphi = (Table(iir,iit,iip+1)-Table(iir,iit,iip-1))/(2*dp)
                
                ! kernel w.r.t. slowness s
                Ks(iir,iit,iip) = TableADJ(iir,iit,iip)*fun(iir,iit,iip)**2
                
                ! kernel w.r.t. anisotropy xi
                if (sqrt(xi(iir,iit,iip)**2+eta(iir,iit,iip)**2) == 0) then
                    Kxi(iir,iit,iip) = TableADJ(iir,iit,iip) &
                      & *(Ttheta**2/(rr(iir)**2)-Tphi**2/(rr(iir)**2*cos(tt(iit))**2))
                    Keta(iir,iit,iip) = TableADJ(iir,iit,iip)*(-2*Ttheta*Tphi/(rr(iir)**2*cos(tt(iit))))
                else
                    Kxi(iir,iit,iip) = TableADJ(iir,iit,iip)*((-gamma*xi(iir,iit,iip)/sqrt(xi(iir,iit,iip)**2 &
                      & +eta(iir,iit,iip)**2))*Tr**2+Ttheta**2/(rr(iir)**2)-Tphi**2/(rr(iir)**2*cos(tt(iit))**2))
                    Keta(iir,iit,iip) = TableADJ(iir,iit,iip)*((-gamma*eta(iir,iit,iip)/sqrt(xi(iir,iit,iip)**2 &
                      & +eta(iir,iit,iip)**2))*Tr**2-2*Ttheta*Tphi/(rr(iir)**2*cos(tt(iit))))
                end if
                ! kernel w.r.t. anisotropy eta
                                                                 
            end do
        end do
    end do

    ! boundary
    do iir=1,nr
        do iit=1,nt
            Ks(iir,iit,1) = 0
            Ks(iir,iit,np) = 0
            Kxi(iir,iit,1) = 0
            Kxi(iir,iit,np) = 0
            Keta(iir,iit,1) = 0
            Keta(iir,iit,np) = 0
        end do
    end do
    do iir=1,nr
        do iip=1,np
            Ks(iir,1,iip) = 0
            Ks(iir,nt,iip) = 0
            Kxi(iir,1,iip) = 0
            Kxi(iir,nt,iip) = 0
            Keta(iir,1,iip) = 0
            Keta(iir,nt,iip) = 0
        end do
    end do
    do iit=1,nt
        do iip=1,np
            Ks(1,iit,iip) = 0
            Ks(nr,iit,iip) = 0
            Kxi(1,iit,iip) = 0
            Kxi(nr,iit,iip) = 0
            Keta(1,iit,iip) = 0
            Keta(nr,iit,iip) = 0
        end do
    end do

end subroutine

subroutine Parameter_Update_Multigrid(rr,tt,pp,nr,nt,np,all_Kernel,nk,invr,invt,invp,ninvr,ninvt,ninvp,ngrid,stepsize,update_value)
    ! forward modeling grid (fine grid)
    integer :: nr,nt,np,nk  ! number of kernels for all parameters
    double precision :: rr(nr),tt(nt),pp(np),dr,dt,dp
    double precision :: update_value(nr,nt,np,nk),para(nr,nt,np,nk),all_Kernel(nr,nt,np,nk)
    
    ! inversion grid (coarse grid)
    integer :: ninvr,ninvt,ninvp,ngrid   ! we have "ngrid" sets of inversion grid 
    double precision :: invr(ngrid,ninvr),invt(ngrid,ninvt),invp(ngrid,ninvp),dinvr,dinvt,dinvp
    double precision :: invkernel(ninvr,ninvt,ninvp,nk)
    double precision :: stepsize
    

    integer :: iir,iit,iip,iik
    integer :: idr,irt,idp,igrid
    double precision :: r1,r2,r3,Linf,pert

    dr=rr(2)-rr(1); dt=tt(2)-tt(1); dp=pp(2)-pp(1)

    
    

    ! 初始化 initiate kernel
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                do iik=1,nk
                    para(iir,iit,iip,iik) = 0
                    update_value(iir,iit,iip,iik) = 0
                end do
            end do
        end do
    end do
    


    
    do igrid=1,ngrid
        dinvr=invr(igrid,2)-invr(igrid,1); dinvt=invt(igrid,2)-invt(igrid,1); dinvp=invp(igrid,2)-invp(igrid,1)

        ! 初始化 initiate kernel
        do iir=1,ninvr
            do iit=1,ninvt
                do iip=1,ninvp
                    do iik=1,nk
                        invkernel(iir,iit,iip,iik) = 0
                    end do
                end do
            end do
        end do
        
        ! smooth the kernel
        do iir=1,nr
            idr = floor((rr(iir)-invr(igrid,1))/dinvr)+1       ! invr(idr) <= rr(iir) < invr(idr+1)
            if (idr<0 .or. idr>ninvr) then
                cycle   ! 如果与反演网格点无关，考察下一个点 if r is out of inversion grid invr, turn to next point
            end if
            r1 = (rr(iir)-invr(igrid,1)-(idr-1)*dinvr)/dinvr   ! r1=0 -> r=invr(idr); r1=1 -> r=invr(idr+1)

            do iit=1,nt
                idt = floor((tt(iit)-invt(igrid,1))/dinvt)+1       
                if (idt<0 .or. idt>ninvt) then
                    cycle   ! 如果与反演网格点无关，考察下一个点 
                end if
                r2 = (tt(iit)-invt(igrid,1)-(idt-1)*dinvt)/dinvt   
                
                do iip=1,np
                    idp = floor((pp(iip)-invp(igrid,1))/dinvp)+1       
                    if (idp<0 .or. idp>ninvp) then
                        cycle   ! 如果与反演网格点无关，考察下一个点 
                    end if
                    r3 = (pp(iip)-invp(igrid,1)-(idp-1)*dinvp)/dinvp   
                    
                    if (r1<0 .or. r1>1 .or. r2<0 .or. r2>1 .or. r3<0 .or. r3>1) then
                        print *, 'error ratio'
                        pause
                    end if
                    
                    do iik=1,nk
                        if (idr>=1 .and. idr<=ninvr .and. idt>=1 .and. idt<=ninvt .and.idp>=1 .and. idp<=ninvp) then
                            invkernel(idr,idt,idp,iik) = invkernel(idr,idt,idp,iik) &
                                                    & + (1-r1)*(1-r2)*(1-r3)*all_Kernel(iir,iit,iip,iik)
                        end if
                        if (idr+1>=1 .and. idr+1<=ninvr .and. idt>=1 .and. idt<=ninvt .and.idp>=1 .and. idp<=ninvp) then
                            invkernel(idr+1,idt,idp,iik) = invkernel(idr+1,idt,idp,iik) &
                                                    & + r1*(1-r2)*(1-r3)*all_Kernel(iir,iit,iip,iik)
                        end if
                        if (idr>=1 .and. idr<=ninvr .and. idt+1>=1 .and. idt+1<=ninvt .and.idp>=1 .and. idp<=ninvp) then
                            invkernel(idr,idt+1,idp,iik) = invkernel(idr,idt+1,idp,iik) &
                                                    & + (1-r1)*r2*(1-r3)*all_Kernel(iir,iit,iip,iik)
                        end if
                        if (idr+1>=1 .and. idr+1<=ninvr .and. idt+1>=1 .and. idt+1<=ninvt .and.idp>=1 .and. idp<=ninvp) then
                            invkernel(idr+1,idt+1,idp,iik) = invkernel(idr+1,idt+1,idp,iik) &
                                                    & + r1*r2*(1-r3)*all_Kernel(iir,iit,iip,iik)
                        end if
                        if (idr>=1 .and. idr<=ninvr .and. idt>=1 .and. idt<=ninvt .and.idp+1>=1 .and. idp+1<=ninvp) then
                            invkernel(idr,idt,idp+1,iik) = invkernel(idr,idt,idp+1,iik) &
                                                    & + (1-r1)*(1-r2)*r3*all_Kernel(iir,iit,iip,iik)
                        end if
                        if (idr+1>=1 .and. idr+1<=ninvr .and. idt>=1 .and. idt<=ninvt .and.idp+1>=1 .and. idp+1<=ninvp) then
                            invkernel(idr+1,idt,idp+1,iik) = invkernel(idr+1,idt,idp+1,iik) &
                                                    & + r1*(1-r2)*r3*all_Kernel(iir,iit,iip,iik)
                        end if
                        if (idr>=1 .and. idr<=ninvr .and. idt+1>=1 .and. idt+1<=ninvt .and.idp+1>=1 .and. idp+1<=ninvp) then
                            invkernel(idr,idt+1,idp+1,iik) = invkernel(idr,idt+1,idp+1,iik) &
                                                    & + (1-r1)*r2*r3*all_Kernel(iir,iit,iip,iik)
                        end if
                        if (idr+1>=1 .and. idr+1<=ninvr .and. idt+1>=1 .and. idt+1<=ninvt .and.idp+1>=1 .and. idp+1<=ninvp) then
                            invkernel(idr+1,idt+1,idp+1,iik) = invkernel(idr+1,idt+1,idp+1,iik) &
                                                    & + r1*r2*r3*all_Kernel(iir,iit,iip,iik)
                        end if

                    
                    end do

                end do
            end do
        end do
        
        ! build update value
        do iir=1,nr
            idr = floor((rr(iir)-invr(igrid,1))/dinvr)+1       ! invr(idr) <= rr(iir) < invr(idr+1)
            if (idr<0 .or. idr>ninvr) then
                cycle   ! 如果与反演网格点无关，考察下一个点 if r is out of inversion grid invr, turn to next point
            end if
            r1 = (rr(iir)-invr(igrid,1)-(idr-1)*dinvr)/dinvr   ! r1=0 -> r=invr(idr); r1=1 -> r=invr(idr+1)

            do iit=1,nt
                idt = floor((tt(iit)-invt(igrid,1))/dinvt)+1       
                if (idt<0 .or. idt>ninvt) then
                    cycle   ! 如果与反演网格点无关，考察下一个点 
                end if
                r2 = (tt(iit)-invt(igrid,1)-(idt-1)*dinvt)/dinvt   
        
                do iip=1,np
                    idp = floor((pp(iip)-invp(igrid,1))/dinvp)+1       
                    if (idp<0 .or. idp>ninvp) then
                        cycle   ! 如果与反演网格点无关，考察下一个点 
                    end if
                    r3 = (pp(iip)-invp(igrid,1)-(idp-1)*dinvp)/dinvp   
                    
                    if (r1<0 .or. r1>1 .or. r2<0 .or. r2>1 .or. r3<0 .or. r3>1) then
                        print *, 'error ratio'
                        pause
                    end if
                    
                    
                    do iik=1,nk
                        pert = 0.0
                        if (idr>=1 .and. idr<=ninvr .and. idt>=1 .and. idt<=ninvt .and.idp>=1 .and. idp<=ninvp) then
                            pert = pert + invkernel(idr,idt,idp,iik)*(1-r1)*(1-r2)*(1-r3) 
                        end if
                        if (idr+1>=1 .and. idr+1<=ninvr .and. idt>=1 .and. idt<=ninvt .and.idp>=1 .and. idp<=ninvp) then
                            pert = pert + invkernel(idr+1,idt,idp,iik)*r1*(1-r2)*(1-r3) 
                        end if
                        if (idr>=1 .and. idr<=ninvr .and. idt+1>=1 .and. idt+1<=ninvt .and.idp>=1 .and. idp<=ninvp) then
                            pert = pert + invkernel(idr,idt+1,idp,iik)*(1-r1)*r2*(1-r3) 
                        end if
                        if (idr+1>=1 .and. idr+1<=ninvr .and. idt+1>=1 .and. idt+1<=ninvt .and.idp>=1 .and. idp<=ninvp) then
                            pert = pert + invkernel(idr+1,idt+1,idp,iik)*r1*r2*(1-r3) 
                        end if
                        if (idr>=1 .and. idr<=ninvr .and. idt>=1 .and. idt<=ninvt .and.idp+1>=1 .and. idp+1<=ninvp) then
                            pert = pert + invkernel(idr,idt,idp+1,iik)*(1-r1)*(1-r2)*r3 
                        end if
                        if (idr+1>=1 .and. idr+1<=ninvr .and. idt>=1 .and. idt<=ninvt .and.idp+1>=1 .and. idp+1<=ninvp) then
                            pert = pert + invkernel(idr+1,idt,idp+1,iik)*r1*(1-r2)*r3 
                        end if
                        if (idr>=1 .and. idr<=ninvr .and. idt+1>=1 .and. idt+1<=ninvt .and.idp+1>=1 .and. idp+1<=ninvp) then
                            pert = pert + invkernel(idr,idt+1,idp+1,iik)*(1-r1)*r2*r3 
                        end if
                        if (idr+1>=1 .and. idr+1<=ninvr .and. idt+1>=1 .and. idt+1<=ninvt .and.idp+1>=1 .and. idp+1<=ninvp) then
                            pert = pert + invkernel(idr+1,idt+1,idp+1,iik)*r1*r2*r3 
                        end if
                        para(iir,iit,iip,iik) = para(iir,iit,iip,iik)+pert
                    end do
                    
                    
                end do
            end do
        end do


    end do

    
    ! rescale
    Linf = 0.0
    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                do iik=1,nk
                    if (Linf < abs(para(iir,iit,iip,iik))) then
                        Linf = abs(para(iir,iit,iip,iik))
                    end if
                end do
            end do
        end do
    end do

    do iir=1,nr
        do iit=1,nt
            do iip=1,np
                do iik=1,nk
                    update_value(iir,iit,iip,iik) = para(iir,iit,iip,iik)/Linf*stepsize
                end do
            end do
        end do
    end do

end subroutine

subroutine Kernel_Mask(rr,tt,pp,nr,nt,np,kernel,rsou,tsou,psou,radius)
    integer :: nr,nt,np
    double precision :: rr(nr),tt(nt),pp(np),rsou,tsou,psou,radius
    double precision :: kernel(nr,nt,np),dis
    integer :: i,j,k

    do i=1,nr
        do j=1,nt
            do k=1,np
                dis = sqrt((rr(i)-rsou)**2+(rsou*(tt(j)-tsou))**2+(rsou*cos(tsou)*(pp(k)-psou))**2)
                if ( dis < radius) then
                    kernel(i,j,k) = kernel(i,j,k) * 0
                end if
            end do
        end do
    end do
end subroutine

subroutine Kernel_Mask_new(rr,tt,pp,nr,nt,np,kernel,rsou,tsou,psou)
    integer :: nr,nt,np
    double precision :: rr(nr),tt(nt),pp(np),rsou,tsou,psou,radius
    double precision :: kernel(nr,nt,np),dis
    integer :: i,j,k
    dr = rr(2)-rr(1); dt = tt(2)-tt(1); dp = pp(2)-pp(1); 
    do i=1,nr
        do j=1,nt
            do k=1,np
                if ( abs(rr(i)-rsou)<dr .and. abs(tt(j)-tsou)<dt .and. abs(pp(k)-psou)<dp) then
                    kernel(i,j,k) = kernel(i,j,k) * 0
                end if
            end do
        end do
    end do
end subroutine

subroutine Adjoint_Source_Dt(Ttime,TtimeT,nrec,sourceADJ,obj)
    integer :: nrec
    double precision :: Ttime(nrec),TtimeT(nrec),sourceADJ(nrec),obj
    integer :: i

    sourceADJ = Ttime-TtimeT
    obj = 0
    do i=1,nrec
        obj = obj + (Ttime(i)-TtimeT(i))**2
    end do

end subroutine

