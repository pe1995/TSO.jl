abstract type RadiativeTransferSolver end

struct ShortCharacteristics{T<:AbstractFloat, I<:Integer} <:RadiativeTransferSolver
    z      ::F
    ntau   ::I
    nbins  ::I
    source ::T
    error  ::T
    dum    ::T
    p      ::T
    sp1    ::T
    sp2    ::T
    sp3    ::T
    an     ::T
    ad     ::T
    bd     ::T
    fact   ::T
    dso    ::T
    dum2   ::T
    a      ::T

    xj     ::T
    x      ::T
    s      ::T
    bplan  ::T
    xh     ::T
    xk     ::T
    ifadim ::I
    kfadim ::I

    pi_s ::Float32
    pi_d ::Float64
    one  ::Float64
    zero ::Float64

    mmy   ::I
    xmy   ::Float32
    xmy2  ::Float32
    wmy   ::Float32
end

function reset!(solver::ShortCharacteristics, κ, S)
    for f in fieldnames(ShortCharacteristics)
        if typeof(f) <: AbstractArray
            getfield(solver, f) .= 0.0
        end
    end

    solver.x     .= κ
    solver.bplan .= S
end

function ShortCharacteristics(T, ntau, nbins)
    z       = zeros(T, ntau)
    source  = zeros(T, ntau)
    error   = zeros(T, ntau)
    dum     = zeros(T, ntau, 3)
    p       = zeros(T, ntau)
    sp1     = zeros(T, ntau, 6)
    sp2     = zeros(T, ntau, 6)
    sp3     = zeros(T, ntau, 6)
    an      = zeros(T, ntau, 6)
    ad      = zeros(T, ntau)
    bd      = zeros(T, ntau)
    fact    = zeros(T, ntau)
    dso     = zeros(T, ntau)
    dum2    = zeros(T, 3, 6)
    a       = zeros(T, 7)
    xj      = zeros(T, ntau)
    x       = zeros(T, ntau)
    s       = zeros(T, ntau)
    bplan   = zeros(T, ntau)
    xh      = zeros(T, ntau)
    xk      = zeros(T, ntau)
    ifadim  = 27000
    kfadim  = 129000
    pi_s    = Base.convert(Float32, 3.141592653589)
    pi_d    = 3.141592653589
    one     = 1.0
    zero    = 0.0
    one_sp  = Base.convert(Float32, 1.0)
    zero_sp = Base.convert(Float32, 0.0)
    mmy     = 0
    xmy     = zeros(Float32, 6)
    xmy2    = zeros(Float32, 6)
    wmy     = zeros(Float32, 6)
    gausi!(mmy, zero_sp, one_sp, wmy, xmy)

    ShortCharacteristics(
        z     ,
        ntau  ,
        nbins ,
        source,
        error ,
        dum   ,
        p     ,
        sp1   ,
        sp2   ,
        sp3   ,
        an    ,
        ad    ,
        bd    ,
        fact  ,
        dso   ,
        dum2  ,
        a     ,
        xj    ,
        x     ,
        s     ,
        bplan ,
        xh    ,
        xk    ,
        ifadim,
        kfadim,
        pi_s  ,
        pi_d  ,
        one   ,
        zero  ,
        mmy   ,
        xmy   ,
        xmy2  ,
        wmy   
        )
end



## Solving

function radiative_transfer(solver, model, eos, opacities)

    # model to perform radiative transfer on
    z, lnρ, lnE = model[:, 1], log.(model[:, 3]), log.(model[:, 2]) 
    T = eltype(z)

    # lookup the optical depth and the source function for this model
    κ  = zeros(T, length(lnρ))
    S  = zeros(T, length(lnρ))
    Iν = zeros(T, size(κ, 2))
    for i in axes(κ, 2)
        κ .= lookup(eos, opacities, :κ,   lnρ, lnE, i)
        S .= lookup(eos, opacities, :src, lnρ, lnE, i)
    
        # solve the radiative transfer
        reset!(solver, κ, S)
        Iν[i] = solve(solver)
    end
end

initiate!(solver::ShortCharacteristics) = begin
    solver.fact   .= 1.0
    solver.dso    .= 0.0
    solver.xj     .= 0.0
    solver.xk     .= 0.0
    solver.error  .= 0.0
    solver.source .= 0.0
end

"""
traneq solves the transfer equation including continuum scattering.
features:
    1. cannons perturbation technique is used on the angular quadrature.
    the basic idea in this technique is to replace the inversion of
    a complicated (mmu order) operator with the inversion of a simple
    operator (one point=eddington approximation), plus iteration on
    the error.

    2. aitken extrapolation accellerates the convergence.

    3. a trick due to robert stein (priv. comm., 1979) is used to
    eliminate the need for double precision storage of the matrix
    elements. the idea is to store the (small) sum of the three
    matrix elements on a row, instead of the (large) diagonal ele-
    ment.

    4. the solution is a cubic spline, rather than a piece-wise
    quadratic function. this is accomplished with the correction
    terms ad and bd in subroutine tranfr.

    5. the scattering is treated as dipole scattering instead of the normally
    used isotropic approximation. this can be done very simply in the
    iterating cannon scheme.

    6. a boundary condition which includes an estimated infalling
    radiation makes the solution good also for values of x+s
    large compared with 1./tau(1). a logarithmic tau-scale
    should be used.

    this version of traneq is compatible with previous traneq's.
    79.06.21 *nord*
"""
function solve(solver::ShortCharacteristics)

    #=  use vartype
    use stateparam
    use ctran
    common /space2/source(mtau),error(mtau),dum(mtau,3),p(mtau)
    & ,sp1(mtau,6),sp2(mtau,6),sp3(mtau,6),an(mtau),ad(mtau),bd(mtau)
    & ,fact(mtau),dso(mtau),dum2(3,6)
    real(sp) :: source,error,dum,p,sp1,sp2,sp3,an,ad,bd,fact,dso,dum2
    real(sp) :: a(7),ds=#

    T     = eltype(solver.x)
    ntau  = solver.ntau
   

    initiate!(solver)


    # calculate the matrix elements
    #call tranfr
    #call transc

    #=
    # iteration loop
    itmax=7
    do 110 it=1,itmax
    110   a(it)=0.
        do 140 it=1,itmax
        itm=it
    c
    c solve the continuum scattering problem in the eddington approximation
        call scattr
        do 120 k=1,jtau
        xj(k)=xj(k)+p(k)
        xk(k)=xk(k)+.333333*p(k)
    c
    c aitken extrapolation used for convergence accelleration
        ds=error(k)+p(k)*s(k)/(x(k)+s(k))
        if(dso(k).gt.0.000001) fact(k)=min(1.25,max(.8,fact(k)-ds/dso(k)))
        ds=ds/fact(k)
        if(it.ge.2) dso(k)=ds
    120   source(k)=source(k)+ds

    c
    c solve the transfer equation with given source function
        call formal
    c
    c check error in source function
        do k=1,jtau
        a(it)=max(a(it),abs(error(k)/source(k)))
        enddo
    c
    c end of iteration loop
        if(a(it).lt.0.0001) go to 141
    140   continue
        write(iwrit,50)(a(it),it=1,itm)
    50    format(' traneq: max error =',12f10.7)
    141   continue
    c
        return
    c---- debug subchk
        end
    =#
    end



"""
formal solves the transfer equation with given source function 'source'.
'error' is the resulting error in the definition of the continuum
scattering source function. transfr calculates the matrix elements
of the problem. flux and intensities at tau=0 are returned in /csurf/.
79.06.21 *nord*
"""
function tranfr(solver)
     
    itran = 0
    c  = zeros(Float32, 6)
    t  = zeros(Float32, 6)
    ex = zeros(Float32, 6)
    y1 = zeros(Float32, 6)
    hsurf ::Float32 = 0.0
    
    a ::Float32
    b ::Float32
    dtaua ::Float32
    dtaub ::Float32
    dtauc ::Float32
    s0 ::Float32
    r1 ::Float32
    p0 ::Float32

    jtau  = solver.ntau
    jtau1 = jtau-1
    jtau2 = jtau-2
    
    for i=1:solver.mmu
        dtaub = 0.5*(solver.x[2] + solver.x[1])*( (solver.z[2]-solver.z[1]) - (solver.z[2]-solver.z[1]))/xmu(i)
        a=1.0/dtaub
        b=a^2
        solver.sp2[1,i] = 1. +2. *a
        solver.sp3[1,i] = -2. *b
        # let p be the even part of the intensity, then
        #
        #         p(2)= p(1) + d*p'(1) + .5*d2*p''(1)
        # or      p(2)= p(1) + d*(p(1)-i(1,-mu)) + .5*d2*(p(1)-s(1)) .
        # where   i(1,-mu) = s(1)*(1.-exp(-t))
        #
        # the difference as compared to the usual second order boundary condition
        # is the additional term   i(1,-mu)=s(1)*(1.-exp(-t)). thus the coefficient
        # for s(1) in the first equation should be changed as follows
        #         s(1)=s(1)*(1.+c*(1.-exp(-t))
        # where   c=2./d
        # *nord* 751009
        c[i]  = 2. *a
        t[i]  = (solver.z[2]-solver.z[1])*(solver.x[2] + solver.x[1])/xmu(i)
        ex[i] = t[i]*(1. -0.5 *t[i] *(1. -.3333 *t[i]))
        if (t[i] > 0.1) 
            ex[i] = 1. -exp(-t[i])
        end

     
            do 100 k=2,jtau1
            dtaua=dtaub
            dtaub=.5*(x(k)+s(k)+x(k+1)+s(k+1))*(tau(k+1)-tau(k))/xmu(i)
            dtauc=.5*(dtaua+dtaub)
            if(itran.eq.0)then
            ad(k)=0.
            bd(k)=0.
            else
            ad(k)=.166667*dtaua/dtauc
            bd(k)=.166667*dtaub/dtauc
            endif
            sp1(k,i)=-1./(dtaua*dtauc)+ad(k)
            sp2(k,i)=1.
        100   sp3(k,i)=-1./(dtaub*dtauc)+bd(k)
        c
        c k=jtau
            sp2(jtau,i)=1.
        c
        c end of mu loop
        110   continue
c
c eliminate subdiagonal, save factors in sp1
      do 121 i=1,mmu
      do 120 k=1,jtau2
      sp1(k,i)=-sp1(k+1,i)/(sp2(k,i)-sp3(k,i))
      sp2(k+1,i)=sp2(k+1,i)+sp1(k,i)*sp2(k,i)
120   sp2(k,i)=sp2(k,i)-sp3(k,i)
121   sp2(jtau-1,i)=sp2(jtau-1,i)-sp3(jtau-1,i)
c
      return
c
      entry formal
c
c zeroset
      do 130 k=1,jtau
      an(k)=(3.*xk(k)-xj(k))/8.*s(k)/(x(k)+s(k))
      xk(k)=0.
130   xj(k)=0.
c
c mu loop
      xh(1)=0.
      hsurf=0.
      do 170 i=1,mmu
c
c initiate approximative source function
      p(1)=source(1)+an(1)*(3.*xmu2(i)-1.)
c note the anisotropic scattering correction
      s0=p(1)
      p(1)=p(1)*(1.+c(i)*ex(i))
      do 140 k=2,jtau1
140   p(k)=(1.-ad(k)-bd(k))*(source(k)+an(k)*(3.*xmu2(i)-1.))
     & +ad(k)*(source(k-1)+an(k-1)*(3.*xmu2(i)-1.))
     & +bd(k)*(source(k+1)+an(k+1)*(3.*xmu2(i)-1.))
      p(jtau)=source(jtau)
c
c accumulate right hand side
      do 150 k=1,jtau2
150   p(k+1)=p(k+1)+sp1(k,i)*p(k)
c
c backsubstitute
      do 160 k=1,jtau1
      p(jtau-k)=(p(jtau-k)/sp2(jtau-k,i)
     .     -(sp3(jtau-k,i)/sp2(jtau-k,i))*p(jtau-k+1))
      xk(jtau-k)=xk(jtau-k)+h(i)*p(jtau-k)*xmu2(i)
160   xj(jtau-k)=xj(jtau-k)+h(i)*p(jtau-k)
c
c end of mu loop
      xk(jtau)=xk(jtau)+h(i)*p(jtau)*xmu2(i)
      r1=p(1)-s0*ex(i)
      xh(1)=xh(1)+h(i)*xmu(i)*r1
      p0=p(1)*(1.-ex(i))+.5*s0*ex(i)**2
      hsurf=hsurf+h(i)*xmu(i)*p0
      y1(i)=2.*p0
c hsurf and y1(6) are the flux and intensities at the surface
170   continue
      xj(jtau)=p(jtau)
c
c 'xj' is the new mean intensity
      do 180 k=1,jtau
180   error(k)=(x(k)*bplan(k)+s(k)*xj(k))/(x(k)+s(k))-source(k)
c
c flux and second moment
      do 190 k=2,jtau
190   xh(k)=2.*(xk(k)-xk(k-1))/(x(k)+s(k)+x(k-1)+s(k-1))/
     /(tau(k)-tau(k-1))
c
      return
c---- debug subchk
      end

function gausi!(k,a,b,ai,xmyi)
    ap    = zeros(29)
    xmyp  = zeros(29)
    indov = zeros(Int, 9)

    ap .= [  1.0, 0.55555555555555,.88888888888888,.347854845137,
             0.65214515486254 ,0.23692688505618, 0.47862867049936,
             0.56888888888888 ,0.17132449237917, 0.36076157304813,
             0.46791393457269 ,0.12948496616887, 0.27970539148927,
             0.38183005050511 ,0.41795918367346, 0.10122853629037,
             0.22238103445337 ,0.31370664587788, 0.36268378337836,
             0.08127438836157 ,0.18064816069485, 0.26061069640293,
             0.31234707704000 ,0.33023935500126, 0.06667134430868,
             0.14945134915058 ,0.21908636251598, 0.26926671930999,
             0.29552422471475 ]
    xmyp .=  [
             0.57735026918962, 0.77459666924148, 0.0,0.86113631159405,
             0.33998104358485, 0.90617984593866, 0.53846931010568, 0.0,
             0.93246951420315, 0.66120938646626, 0.23861918608319,
             0.94910791234275, 0.74153118559939, 0.40584515137739, 0.0,
             0.96028985649753, 0.79666647741362, 0.52553240991632,
             0.18343464249565, 0.96816023950762, 0.83603110732663,
             0.61337143270059, 0.32425342340380, 0.0,0.97390652851717,
             0.86506336668898, 0.67940956829902, 0.43339539412924,
             0.14887433898163]

    indov = [1,3,5,8,11,15,19,24,29]

    if (k == 1)
        xmyi[1]= (b+a)*0.5
        ai[1]  = b-a
    else
        kud = 0
        flk = float(k)/2.
        k2  = k/2
        fk  = float(k2)
        if (abs(flk-fk)-1.e-7 < 0.0)
            nothing
        else
            k2=k2+1
            kud=1
        end

        ioev=indov[k-1]
        ined=ioev-k2
        for i=1:k2
            ip      = ined+i
            xmyi[i] = -xmyp[ip]*(b-a)*0.5+(b+a)*0.5
            ai[i]   = (b-a)*0.5*ap[ip]
        end
        k2=k2+1
        for i=k2,k
            ip=ioev+k2-i
            if (kud<=0)
                nothing
            else
                ip=ip-1
            end
            xmyi[i] = xmyp[ip]*(b-a)*0.5+(b+a)*0.5
            ai[i]   = (b-a)*0.5*ap[ip]
        end
    end
end