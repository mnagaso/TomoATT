!  3-D third order LF scheme (point source, sphere, 2nd order accuracy, parallel) (verified by ega1)
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
    double precision,parameter :: tol = (10.0)**(-4),eps=10.0**(-12)
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
                        print *, rr(iir),t0*180/3.1415927,iit,p0*180/3.1415927,iip
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
                write(*,*) 'id',myid,'iter ',iter,', T is steadt at L1 dif= ',L1_dif, ', Linf dif= ',Linf_dif
                isexit = .true.
            else
            !    write(*,*) 'id',myid,'iter ',iter,', T is changing, L1 dif = ', L1_dif,'L inf dif = ', Linf_dif
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

! multiplicative factors are more robust
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
    double precision,parameter :: tol = (10.0)**(-4),eps=10.0**(-12)
    integer,parameter :: MaxIter=500
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

    !open(10,file='matlab/ega2/160_160_160/timetable_mul')
    !do iir=3,nr-2
    !    do iit=3,nt-2
    !        do iip=3,np-2
    !            read(10,*) tmp
    !            tau(iir,iit,iip) = tmp/T0v(iir,iit,iip)
    !        end do
    !    end do
    !end do
    !close(10)

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
    write(*,*) 'T0 L1 error is ',L1_err,', T0 Linf error is', Linf_err




    ! 正式迭代，更新tau
    do iter =1,MaxIter
        tau_old = tau
        L1_dif=0; Linf_dif=0;L1_err=0;Linf_err=0;
        do rdirec = -1,1,2
            do tdirec = -1,1,2
                do pdirec = -1,1,2

                    !x: nr-1 <-> 2, y: nt-1 <-> 2, z: np-1 <-> 2
                    !sigr = 0; sigt=0; sigp=0
                    !do iir=1,nr
                    !    do iit=1,nt
                    !        do iip=1,np
                    !            sigr = max(sigr,2*sqrt(a(iir,iit,iip))*T0v(iir,iit,iip))
                    !            sigt = max(sigt,2*sqrt(b(iir,iit,iip))*T0v(iir,iit,iip))
                    !            sigp = max(sigp,2*sqrt(c(iir,iit,iip))*T0v(iir,iit,iip))
                    !        end do
                    !    end do
                    !end do
                    !coe=1.0/((sigr/dr)+(sigt/dt)+(sigp/dp))


                    do iir=nint(0.5+nr/2.0+(nr/2.0-1.5)*rdirec),nint(0.5+nr/2.0+(-nr/2.0+1.5)*rdirec),-rdirec
                        do iit=nint(0.5+nt/2.0+(nt/2.0-1.5)*tdirec),nint(0.5+nt/2.0+(-nt/2.0+1.5)*tdirec),-tdirec
                            do iip=nint(0.5+np/2.0+(np/2.0-1.5)*pdirec),nint(0.5+np/2.0+(-np/2.0+1.5)*pdirec),-pdirec

                                if(ischange(iir,iit,iip)==1) then
                                    sigr = 2*sqrt(a(iir,iit,iip))*T0v(iir,iit,iip)
                                    sigt = 2*sqrt(b(iir,iit,iip))*T0v(iir,iit,iip)
                                    sigp = 2*sqrt(c(iir,iit,iip))*T0v(iir,iit,iip)
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
                            !if (tau(3,iit,iip)>0) then
                                tau(1,iit,iip) = max(2*tau(2,iit,iip)-tau(3,iit,iip),tau(3,iit,iip))
                            !else
                            !    tau(1,iit,iip) = min(2*tau(2,iit,iip)-tau(3,iit,iip),tau(3,iit,iip))
                            !end if
                            !if (tau(nr-2,iit,iip)>0) then
                                tau(nr,iit,iip) = max(2*tau(nr-1,iit,iip)-tau(nr-2,iit,iip),tau(nr-2,iit,iip))
                            !else
                            !    tau(nr,iit,iip) = min(2*tau(nr-1,iit,iip)-tau(nr-2,iit,iip),tau(nr-2,iit,iip))
                            !end if
                        end do
                    end do
                    do iir=1,nr
                        do iip=1,np
                            !if (tau(iir,3,iip)>0) then
                                tau(iir,1,iip) = max(2*tau(iir,2,iip)-tau(iir,3,iip),tau(iir,3,iip))
                            !else
                            !    tau(iir,1,iip) = min(2*tau(iir,2,iip)-tau(iir,3,iip),tau(iir,3,iip))
                            !end if
                            !if (tau(iir,nt-2,iip)>0) then
                                tau(iir,nt,iip) = max(2*tau(iir,nt-1,iip)-tau(iir,nt-2,iip),tau(iir,nt-2,iip))
                            !else
                            !    tau(iir,nt,iip) = min(2*tau(iir,nt-1,iip)-tau(iir,nt-2,iip),tau(iir,nt-2,iip))
                            !end if
                        end do
                    end do
                    do iir=1,nr
                        do iit=1,nt
                            !if (tau(iir,iit,3)>0) then
                                tau(iir,iit,1) = max(2*tau(iir,iit,2)-tau(iir,iit,3),tau(iir,iit,3))
                            !else
                            !    tau(iir,iit,1) = min(2*tau(iir,iit,2)-tau(iir,iit,3),tau(iir,iit,3))
                            !end if
                            !if (tau(iir,iit,np-2)>0) then
                                tau(iir,iit,np) = max(2*tau(iir,iit,np-1)-tau(iir,iit,np-2),tau(iir,iit,np-2))
                            !else
                            !    tau(iir,iit,np) = min(2*tau(iir,iit,np-1)-tau(iir,iit,np-2),tau(iir,iit,np-2))
                            !end if
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

        if (abs(L1_dif)<tol ) then
            write(*,*) 'iter ',iter,', T is steadt'
            exit
        else
            write(*,*) 'iter ',iter,', T is changing, continue ... '
        end if

        if (iter==MaxIter) then
            write(*,*) 'iter ',iter,', max iteration steps'
        end if

        write(*,'(a,f15.7,a,f15.7)') 'L_1(T-u)=',L1_err,'  L_inf(T_LF-u)',Linf_err

    end do

    T=T0v*tau

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
subroutine Linear_Interp_3D(xx,yy,zz,val,nx,ny,nz,x0,y0,z0,v0,n0)
    integer :: nx,ny,nz,idx0,idy0,idz0,n0
    double precision :: xx(nx),yy(ny),zz(nz),val(nx,ny,nz),x0(n0),y0(n0),z0(n0),v0(n0),dx,dy,dz
    double precision :: r1,r2,r3
    integer :: iix,iiy,iiz,i

    dx = xx(2)-xx(1)
    dy = yy(2)-yy(1)
    dz = zz(2)-zz(1)

    do i=1,n0
        idx0 = floor((x0(i)-xx(1))/dx)+1
        idy0 = floor((y0(i)-yy(1))/dy)+1
        idz0 = floor((z0(i)-zz(1))/dz)+1

        if (idx0<1 .or. idx0>nx-1 .or. idy0<1 .or. idy0>ny-1 .or. idz0<1 .or. idz0>nz-1) then
            write(*,*) 'point out of the mesh'
            pause
            return
        end if

        r1 = min(1.0, (x0(i)-xx(idx0))/dx )
        r2 = min(1.0, (y0(i)-yy(idy0))/dy )
        r3 = min(1.0, (z0(i)-zz(idz0))/dz )

        v0(i) = (1-r1)*(1-r2)*(1-r3)*val(idx0,idy0,idz0) + (1-r1)*(1-r2)*r3*val(idx0,idy0,idz0+1) &
            & + (1-r1)*r2*(1-r3)*val(idx0,idy0+1,idz0) + (1-r1)*r2*r3*val(idx0,idy0+1,idz0+1) &
            & + r1*(1-r2)*(1-r3)*val(idx0+1,idy0,idz0) + r1*(1-r2)*r3*val(idx0+1,idy0,idz0+1) &
            & + r1*r2*(1-r3)*val(idx0+1,idy0+1,idz0) + r1*r2*r3*val(idx0+1,idy0+1,idz0+1)
    end do
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
                !tmpr2 = (rr(iir)+rr(iir+1))/2; tmpt2 = (tt(iit)+tt(iit-1))/2;
                tmpr2 = (rr(iir)+rr(iir+1))/2; tmpt2 = (tt(iit)+tt(iit+1))/2;

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
            write(*,*) 'iter ',iter,', TableADJ is steady at ',Linf_dif
            exit
        else
        !    write(*,*) 'iter ',iter,', TableADJ is changing, continue ... '
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

subroutine Adjoint_Source_Ddt(Ttime,TtimeT,nrec,sourceADJ,obj)
    integer :: nrec
    double precision :: Ttime(nrec),TtimeT(nrec),sourceADJ(nrec),obj
    integer :: i,j

    sourceADJ = 0
    obj = 0
    do i=1,nrec-1
        do j=i,nrec
            obj = obj + ((Ttime(i)-Ttime(j))-(TtimeT(i)-TtimeT(j)))**2
            sourceADJ(i) = sourceADJ(i)+((Ttime(i)-Ttime(j))-(TtimeT(i)-TtimeT(j)))
            sourceADJ(j) = sourceADJ(j)-((Ttime(i)-Ttime(j))-(TtimeT(i)-TtimeT(j)))
        end do
    end do

end subroutine

subroutine Epicenter_Distance(lat1,lon1,lat2,lon2,dis)
    double precision :: lat1,lon1,lat2,lon2,dis
    double precision,parameter :: pi=3.1415926535897932384626433, R=6378.137

    if((lat1-lat2)**2+(lon1-lon2)**2<0.0001) then
        dis=0.0
    else
        dis = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))*R
    end if
end subroutine
