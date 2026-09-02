/* CUDA execution plumbing for AMD2.  Coefficients and diffusion are evaluated
 * on the device, and this source is compiled with strict CUDA arithmetic flags.
 */
#include "boundary.h"
#include "cuda_launcher.h"
#include "diff_amd2_kernels.cuh"
#include "diff_amd2.h"
#include "diff_kl_kernels.cuh"
#include "fields.h"
#include "grid.h"
#include "stats.h"
#include "thermo.h"
#include "tools.h"
#include <limits>
#include <type_traits>

namespace
{
    __device__ void atomic_max_nonnegative(double* address,const double value)
    {
        auto* bits=reinterpret_cast<unsigned long long*>(address);
        unsigned long long old=*bits;
        while (__longlong_as_double(old)<value)
        {
            const unsigned long long assumed=old;
            old=atomicCAS(bits,assumed,__double_as_longlong(value));
            if (old==assumed) break;
        }
    }

    __global__ void reset_aggregates_g(unsigned long long* counts,const int ncounts,double* maximum)
    {
        const int n=blockIdx.x*blockDim.x+threadIdx.x;
        if (n<ncounts) counts[n]=0;
        if (n==0) maximum[0]=0.;
    }

    template<typename TF> __device__ double dz_centered(const TF* f,const int n,const int k,const int kk,const double alpha,const TF* dzhi)
    {
        const double lo=(double(f[n])-double(f[n-kk]))*double(dzhi[k]);
        const double hi=(double(f[n+kk])-double(f[n]))*double(dzhi[k+1]);
        const double low=(1.-alpha)*lo, high=alpha*hi;
        return low+high;
    }
    template<typename TF> __device__ double dz_onesided(const TF* f,const int n,const int kk,const TF* z,const int k0,const int k1,const int k2)
    {
        const double x0=z[k0],x1=z[k1],x2=z[k2];
        const double w0=(2.*x0-x1-x2)/((x0-x1)*(x0-x2));
        const double w1=(x0-x2)/((x1-x0)*(x1-x2));
        const double w2=(x0-x1)/((x2-x0)*(x2-x1));
        return w0*double(f[n])+w1*double(f[n+(k1-k0)*kk])+w2*double(f[n+(k2-k0)*kk]);
    }
    template<typename TF> __device__ double dx_centered(const TF* f,const int n,const double dxi)
    { const double a=double(f[n+1])-double(f[n]),b=double(f[n])-double(f[n-1]); return (a+b)*(0.5*dxi); }
    template<typename TF> __device__ double dy_centered(const TF* f,const int n,const int jj,const double dyi)
    { const double a=double(f[n+jj])-double(f[n]),b=double(f[n])-double(f[n-jj]); return (a+b)*(0.5*dyi); }

    template<typename TF> __device__ bool make_A(
        double A[3][3],double Ah[3][3],double& scale,double& Du,
        const TF* u,const TF* v,const TF* w,const TF* z,const TF* zh,const TF* dzi,const TF* dzhi,
        const double dxi,const double dyi,const int n,const int k,const int jj,const int kk,
        const int kstart,const int kend,const bool ustar_bot,const bool ustar_top)
    {
        const double alpha=(double(z[k])-double(zh[k]))/(double(zh[k+1])-double(zh[k]));
        A[0][0]=(double(u[n+1])-double(u[n]))*dxi;
        A[1][1]=(double(v[n+jj])-double(v[n]))*dyi;
        A[2][2]=(double(w[n+kk])-double(w[n]))*double(dzi[k]);
        const double u0=double(u[n+jj])-double(u[n-jj]),u1=double(u[n+1+jj])-double(u[n+1-jj]);
        A[0][1]=(u0+u1)*(0.25*dyi);
        const double v0=double(v[n+1])-double(v[n-1]),v1=double(v[n+jj+1])-double(v[n+jj-1]);
        A[1][0]=(v0+v1)*(0.25*dxi);
        auto vertical=[&](const TF* f,const int q)
        {
            if (k==kstart&&ustar_bot) return dz_onesided(f,q,kk,z,k,k+1,k+2);
            if (k==kend-1&&ustar_top) return dz_onesided(f,q,kk,z,k,k-1,k-2);
            return dz_centered(f,q,k,kk,alpha,dzhi);
        };
        const double uz0=vertical(u,n),uz1=vertical(u,n+1); A[0][2]=0.5*(uz0+uz1);
        const double vz0=vertical(v,n),vz1=vertical(v,n+jj); A[1][2]=0.5*(vz0+vz1);
        const double wx0=dx_centered(w,n,dxi),wx1=dx_centered(w,n+kk,dxi); A[2][0]=(1.-alpha)*wx0+alpha*wx1;
        const double wy0=dy_centered(w,n,jj,dyi),wy1=dy_centered(w,n+kk,jj,dyi); A[2][1]=(1.-alpha)*wy0+alpha*wy1;
        scale=0.; bool ok=isfinite(alpha);
        for(int i=0;i<3;++i)for(int j=0;j<3;++j){ok=ok&&isfinite(A[i][j]);scale=fmax(scale,fabs(A[i][j]));}
        if(!ok||!isfinite(scale)||scale==0.) return false;
        Du=0.;for(int i=0;i<3;++i)for(int j=0;j<3;++j){Ah[i][j]=A[i][j]/scale;Du+=Ah[i][j]*Ah[i][j];}
        return isfinite(Du)&&Du>0.;
    }

    template<typename TF> __global__ void amd_momentum_g(
        TF* evisc,TF* nus,TF* nub,unsigned long long* counts,double* maximum,
        const TF* u,const TF* v,const TF* w,const TF* b,const TF* bh,
        const TF* z,const TF* zh,const TF* dz,const TF* dzi,const TF* dzhi,
        const double dxi,const double dyi,const double cx,const double cy,const double cz,const double cap,const double tfmax,
        const double d0,const double d1,const double d2,const double G0,const double G1,const double G2,
        const bool use_buoy,const bool ustar_bot,const bool ustar_top,
        const int istart,const int iend,const int jstart,const int jend,const int kstart,const int kend,const int jj,const int kk)
    {
        const double direction[3]={d0,d1,d2},background[3]={G0,G1,G2};
        const int nx=iend-istart,ny=jend-jstart,nz=kend-kstart;
        const int cell=blockIdx.x*blockDim.x+threadIdx.x;
        if(cell>=nx*ny*nz)return;
        const int i=istart+cell%nx;
        const int j=jstart+(cell/nx)%ny;
        const int k=kstart+cell/(nx*ny);
            const int n=i+j*jj+k*kk; evisc[n]=TF(0);nus[n]=TF(0);nub[n]=TF(0);
            double A[3][3],Ah[3][3],a,Du; int status=0;
            bool finite_A=true; // distinguish all-zero from invalid below
            const double alpha=(double(z[k])-double(zh[k]))/(double(zh[k+1])-double(zh[k]));
            A[0][0]=(double(u[n+1])-double(u[n]))*dxi;A[1][1]=(double(v[n+jj])-double(v[n]))*dyi;A[2][2]=(double(w[n+kk])-double(w[n]))*double(dzi[k]);
            const double ux0=double(u[n+jj])-double(u[n-jj]),ux1=double(u[n+1+jj])-double(u[n+1-jj]);A[0][1]=(ux0+ux1)*(0.25*dyi);
            const double vx0=double(v[n+1])-double(v[n-1]),vx1=double(v[n+jj+1])-double(v[n+jj-1]);A[1][0]=(vx0+vx1)*(0.25*dxi);
            auto vertical=[&](const TF* f,const int q){if(k==kstart&&ustar_bot)return dz_onesided(f,q,kk,z,k,k+1,k+2);if(k==kend-1&&ustar_top)return dz_onesided(f,q,kk,z,k,k-1,k-2);return dz_centered(f,q,k,kk,alpha,dzhi);};
            const double uz0=vertical(u,n),uz1=vertical(u,n+1);A[0][2]=0.5*(uz0+uz1);const double vz0=vertical(v,n),vz1=vertical(v,n+jj);A[1][2]=0.5*(vz0+vz1);
            const double wx0=dx_centered(w,n,dxi),wx1=dx_centered(w,n+kk,dxi);A[2][0]=(1.-alpha)*wx0+alpha*wx1;const double wy0=dy_centered(w,n,jj,dyi),wy1=dy_centered(w,n+kk,jj,dyi);A[2][1]=(1.-alpha)*wy0+alpha*wy1;
            a=0.;for(int p=0;p<3;++p)for(int q=0;q<3;++q){finite_A=finite_A&&isfinite(A[p][q]);a=fmax(a,fabs(A[p][q]));}
            if(!finite_A||!isfinite(a)){atomicAdd(&counts[0],1ULL);return;} if(a==0.){atomicAdd(&counts[1],1ULL);return;}
            Du=0.;double S[3][3];for(int p=0;p<3;++p)for(int q=0;q<3;++q){Ah[p][q]=A[p][q]/a;Du+=Ah[p][q]*Ah[p][q];}
            if(!isfinite(Du)||!(Du>0.)){atomicAdd(&counts[0],1ULL);return;}for(int p=0;p<3;++p)for(int q=0;q<3;++q)S[p][q]=0.5*(Ah[p][q]+Ah[q][p]);
            const double M[3]={cx/(dxi*dxi),cy/(dyi*dyi),cz*double(dz[k])*double(dz[k])};double shear=0.;for(int q=0;q<3;++q)for(int p=0;p<3;++p)for(int r=0;r<3;++r)shear+=M[q]*Ah[p][q]*Ah[r][q]*S[p][r];const double ns=-a*shear/Du;
            if(!isfinite(shear)||!isfinite(ns)){atomicAdd(&counts[2],1ULL);return;}double nb=0.;bool bok=true,zero=false;
            if(use_buoy){double qv[3]={dx_centered(b,n,dxi)+background[0],dy_centered(b,n,jj,dyi)+(background[1]),(double(bh[n+kk])-double(bh[n]))*double(dzi[k])+background[2]};double c=0.;for(int q=0;q<3;++q){bok=bok&&isfinite(qv[q]);c=fmax(c,fabs(qv[q]));}if(!isfinite(c))bok=false;else if(c==0.)zero=true;else{double sum=0.;for(int q=0;q<3;++q){double projection=0.;for(int p=0;p<3;++p)projection+=direction[p]*Ah[p][q];sum+=M[q]*projection*(qv[q]/c);}nb=(c/a)*sum/Du;bok=bok&&isfinite(sum)&&isfinite(nb);}}
            if(zero)atomicAdd(&counts[10],1ULL);if(!bok){atomicAdd(&counts[3],1ULL);return;}if(fabs(ns)>tfmax||fabs(nb)>tfmax){atomicAdd(&counts[4],1ULL);return;}const double raw=ns+nb;if(!isfinite(raw)){atomicAdd(&counts[5],1ULL);return;}if(raw<=0.){nus[n]=TF(ns);nub[n]=TF(nb);atomicAdd(&counts[6],1ULL);return;}const double stored=(cap>0.&&raw>cap)?cap:raw;if(!isfinite(stored)||stored>tfmax){atomicAdd(&counts[7],1ULL);return;}nus[n]=TF(ns);nub[n]=TF(nb);evisc[n]=TF(stored);status=stored!=raw?9:8;atomicAdd(&counts[status],1ULL);atomic_max_nonnegative(maximum,stored);
    }

    template<typename TF> __global__ void amd_scalar_g(
        TF* ediff,unsigned long long* counts,double* maximum,const TF* phi,const TF* u,const TF* v,const TF* w,
        const TF* z,const TF* zh,const TF* dz,const TF* dzi,const TF* dzhi,const double dxi,const double dyi,
        const double cx,const double cy,const double cz,const double cap,const double tfmax,const bool ustar_bot,const bool ustar_top,
        const int istart,const int iend,const int jstart,const int jend,const int kstart,const int kend,const int jj,const int kk)
    {
        const int nx=iend-istart,ny=jend-jstart,nz=kend-kstart;
        const int cell=blockIdx.x*blockDim.x+threadIdx.x;
        if(cell>=nx*ny*nz)return;
        const int i=istart+cell%nx,j=jstart+(cell/nx)%ny,k=kstart+cell/(nx*ny);
        const int n=i+j*jj+k*kk;ediff[n]=TF(0);double A[3][3],Ah[3][3],a,Du;if(!make_A(A,Ah,a,Du,u,v,w,z,zh,dzi,dzhi,dxi,dyi,n,k,jj,kk,kstart,kend,ustar_bot,ustar_top))return;const double alpha=(double(z[k])-double(zh[k]))/(double(zh[k+1])-double(zh[k]));double grad[3]={dx_centered(phi,n,dxi),dy_centered(phi,n,jj,dyi),dz_centered(phi,n,k,kk,alpha,dzhi)};double g=0.;bool ok=true;for(int q=0;q<3;++q){ok=ok&&isfinite(grad[q]);g=fmax(g,fabs(grad[q]));}int status=1;if(!ok||!isfinite(g)){atomicAdd(&counts[1],1ULL);return;}if(g==0.){atomicAdd(&counts[0],1ULL);return;}double p[3],D=0.;for(int q=0;q<3;++q){p[q]=grad[q]/g;D+=p[q]*p[q];}if(!isfinite(D)||!(D>0.)){atomicAdd(&counts[1],1ULL);return;}const double M[3]={cx/(dxi*dxi),cy/(dyi*dyi),cz*double(dz[k])*double(dz[k])};double sum=0.;for(int q=0;q<3;++q)for(int r=0;r<3;++r)sum+=M[q]*Ah[r][q]*p[q]*p[r];const double raw=-a*sum/D;if(!isfinite(sum)||!isfinite(raw)){atomicAdd(&counts[2],1ULL);return;}if(raw<=0.){atomicAdd(&counts[3],1ULL);return;}const double stored=(cap>0.&&raw>cap)?cap:raw;if(!isfinite(stored)||stored>tfmax){atomicAdd(&counts[4],1ULL);return;}ediff[n]=TF(stored);status=stored!=raw?6:5;atomicAdd(&counts[status],1ULL);atomic_max_nonnegative(maximum,stored);
    }

    template<typename TF> __global__ void reflect_g(TF* f,const int icells,const int jcells,const int kstart,const int kend,const int kk)
    {const int i=blockIdx.x*blockDim.x+threadIdx.x,j=blockIdx.y*blockDim.y+threadIdx.y;if(i<icells&&j<jcells){const int b=i+j*icells+kstart*kk,t=i+j*icells+(kend-1)*kk;f[b-kk]=f[b];f[t+kk]=f[t];}}

    __global__ void timestep_init_g(unsigned long long* invalid,double* maximum)
    {if(!blockIdx.x&&!threadIdx.x){invalid[0]=0;invalid[1]=0;maximum[0]=0.;maximum[1]=0.;maximum[2]=0.;maximum[3]=0.;}}
    template<typename TF> __global__ void timestep_variable_g(
        unsigned long long* invalid,double* maximum,const TF* coeff,const TF molecular,const TF* dzi,
        const double dx2i,const double dy2i,const int istart,const int iend,const int jstart,const int jend,
        const int kstart,const int kend,const int jj,const int kk)
    {if(blockIdx.x||threadIdx.x)return;for(int k=kstart;k<kend;++k)for(int j=jstart;j<jend;++j)for(int i=istart;i<iend;++i){const int n=i+j*jj+k*kk;const double sgs=coeff?double(coeff[n]):0.;const double H=dx2i+dy2i+double(dzi[k])*double(dzi[k]);const double value=(sgs+double(molecular))*H;if(!isfinite(sgs)||sgs<0.||!isfinite(value)||value<0.){invalid[0]=1;invalid[1]=static_cast<unsigned long long>(n);maximum[1]=sgs;maximum[2]=double(molecular);maximum[3]=H;return;}maximum[0]=fmax(maximum[0],value);}}

}

template<typename TF>
void Diff_amd2<TF>::prepare_device(Boundary<TF>&)
{
    compact_counts_g.allocate(Mom_status_count+1);
    compact_max_g.allocate(4);
}

template<typename TF>
void Diff_amd2<TF>::clear_device()
{
    compact_counts_g.free();
    compact_max_g.free();
}

template<typename TF>
void Diff_amd2<TF>::exec_viscosity(Stats<TF>&, Thermo<TF>& th)
{
    auto& gd=grid.get_grid_data();
    const bool use_buoy=swamd_buoy&&th.get_switch()!=Thermo_type::Disabled;
    std::shared_ptr<Field3d<TF>> b,bh;
    if(use_buoy){b=fields.get_tmp_g();bh=fields.get_tmp_g();th.get_thermo_field_g(*b,"b",true);th.get_thermo_field_g(*bh,"b_h",false);}
    const auto geometry=th.get_buoyancy_geometry();
    const bool ustar_bot=boundary.get_momentum_bcbot()==Boundary_type::Ustar_type;
    const bool ustar_top=boundary.get_momentum_bctop()==Boundary_type::Ustar_type;
    constexpr int coefficient_block_size=256;
    const int coefficient_cells=gd.imax*gd.jmax*gd.kmax;
    const int coefficient_blocks=(coefficient_cells+coefficient_block_size-1)/coefficient_block_size;
    reset_aggregates_g<<<1,Mom_status_count+1>>>(compact_counts_g,Mom_status_count+1,compact_max_g);
    amd_momentum_g<TF><<<coefficient_blocks,coefficient_block_size>>>(
        fields.sd.at("evisc")->fld_g,fields.sd.at("nu_s")->fld_g,fields.sd.at("nu_b")->fld_g,
        compact_counts_g,compact_max_g,fields.mp.at("u")->fld_g,fields.mp.at("v")->fld_g,fields.mp.at("w")->fld_g,
        use_buoy?b->fld_g.data():nullptr,use_buoy?bh->fld_g.data():nullptr,gd.z_g,gd.zh_g,gd.dz_g,gd.dzi_g,gd.dzhi_g,
        double(gd.dxi),double(gd.dyi),camd[0],camd[1],camd[2],amd_max,double(std::numeric_limits<TF>::max()),
        double(geometry.force_direction[0]),double(geometry.force_direction[1]),double(geometry.force_direction[2]),
        double(geometry.background_gradient[0]),double(geometry.background_gradient[1]),double(geometry.background_gradient[2]),
        use_buoy,ustar_bot,ustar_top,gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
    cuda_check_error();
    std::vector<unsigned long long> host_counts(Mom_status_count+1);double host_max;
    cudaMemcpy(host_counts.data(),compact_counts_g,(Mom_status_count+1)*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
    cudaMemcpy(&host_max,compact_max_g,sizeof(double),cudaMemcpyDeviceToHost);
    Aggregates pass;for(int n=0;n<Mom_status_count;++n){pass.mom[n]=static_cast<std::uint64_t>(host_counts[n]);pass.cells_evaluated+=pass.mom[n];}pass.zero_buoyancy_gradient=static_cast<std::uint64_t>(host_counts[Mom_status_count]);pass.max_evisc=host_max;

    dim3 block(gd.ithread_block,gd.jthread_block,1),grid2((gd.icells+block.x-1)/block.x,(gd.jcells+block.y-1)/block.y,1);
    auto halos=[&](Field3d<TF>& f){boundary_cyclic.exec_g(f.fld_g);reflect_g<TF><<<grid2,block>>>(f.fld_g,gd.icells,gd.jcells,gd.kstart,gd.kend,gd.ijcells);};
    halos(*fields.sd.at("evisc"));halos(*fields.sd.at("nu_s"));halos(*fields.sd.at("nu_b"));
    for(const auto& item:scalar_coeff)
    {
        reset_aggregates_g<<<1,Scalar_status_count>>>(compact_counts_g,Scalar_status_count,compact_max_g);
        amd_scalar_g<TF><<<coefficient_blocks,coefficient_block_size>>>(fields.sd.at(item.second)->fld_g,compact_counts_g,compact_max_g,fields.sp.at(item.first)->fld_g,
            fields.mp.at("u")->fld_g,fields.mp.at("v")->fld_g,fields.mp.at("w")->fld_g,gd.z_g,gd.zh_g,gd.dz_g,gd.dzi_g,gd.dzhi_g,
            double(gd.dxi),double(gd.dyi),camd[0],camd[1],camd[2],amd_max,double(std::numeric_limits<TF>::max()),ustar_bot,ustar_top,
            gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
        cudaMemcpy(host_counts.data(),compact_counts_g,Scalar_status_count*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
        cudaMemcpy(&host_max,compact_max_g,sizeof(double),cudaMemcpyDeviceToHost);
        auto& c=pass.scalar[item.first];c.resize(Scalar_status_count);for(int n=0;n<Scalar_status_count;++n)c[n]=static_cast<std::uint64_t>(host_counts[n]);pass.max_scalar[item.first]=host_max;
        halos(*fields.sd.at(item.second));
    }
    cuda_check_error();
    if(use_buoy){fields.release_tmp_g(b);fields.release_tmp_g(bh);}
    add_aggregates(active,pass);add_aggregates(cumulative,pass);
}

template<typename TF>
void Diff_amd2<TF>::exec(Stats<TF>& stats)
{
    auto& gd=grid.get_grid_data();
    const Grid_layout layout={gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.istride,gd.jstride,gd.kstride};
    const TF dxidxi=TF(1)/(gd.dx*gd.dx), dyidyi=TF(1)/(gd.dy*gd.dy);
    auto run=[&](auto surface_tag)
    {
        constexpr bool surface=decltype(surface_tag)::value;
        launch_grid_kernel<Diff_les_kernels::diff_uvw_g<TF,surface,true>>(
            layout,fields.mt.at("u")->fld_g.view(),fields.mt.at("v")->fld_g.view(),fields.mt.at("w")->fld_g.view(),
            fields.sd.at("evisc")->fld_g,fields.mp.at("u")->fld_g,fields.mp.at("v")->fld_g,fields.mp.at("w")->fld_g,
            fields.mp.at("u")->flux_bot_g,fields.mp.at("u")->flux_top_g,fields.mp.at("v")->flux_bot_g,fields.mp.at("v")->flux_top_g,
            gd.dzi_g,gd.dzhi_g,gd.dxi,gd.dyi,fields.rhoref_g,fields.rhorefh_g,fields.rhorefi_g,fields.rhorefhi_g,fields.visc);
        for (const auto& item:fields.st)
        {
            auto& p=*fields.sp.at(item.first);
            if (swamd_scalar)
                launch_grid_kernel<Diff_les_kernels::diff_c_g<TF,surface,true>>(
                    layout,item.second->fld_g.view(),p.fld_g,fields.sd.at(scalar_coeff.at(item.first))->fld_g,
                    p.flux_bot_g,p.flux_top_g,gd.dzi_g,gd.dzhi_g,dxidxi,dyidyi,
                    fields.rhorefi_g,fields.rhorefh_g,TF(1),p.visc);
            else
            {
                launch_grid_kernel<Diff_amd2_kernels::diff_c_molecular_g<TF,surface>>(
                    layout,item.second->fld_g.view(),p.fld_g,p.flux_bot_g,p.flux_top_g,
                    gd.dzi_g,gd.dzhi_g,fields.rhorefi_g,fields.rhorefh_g,dxidxi,dyidyi,p.visc);
            }
        }
    };
    if (boundary.get_switch()=="default") run(std::false_type{}); else run(std::true_type{});
    cuda_check_error(); cudaDeviceSynchronize();
    for (const auto& name:{"u","v","w"}) stats.calc_tend(*fields.mt.at(name),"diff");
    for (const auto& item:fields.st) stats.calc_tend(*item.second,"diff");
}

template<typename TF>
double Diff_amd2<TF>::get_dn(const double dt)
{
    auto& gd=grid.get_grid_data();
    timestep_init_g<<<1,1>>>(compact_counts_g,compact_max_g);
    auto check=[&](const std::string& name,const TF* coeff,const TF molecular)
    {
        timestep_variable_g<TF><<<1,1>>>(compact_counts_g,compact_max_g,coeff,molecular,gd.dzi_g,
            1./(double(gd.dx)*double(gd.dx)),1./(double(gd.dy)*double(gd.dy)),gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
        unsigned long long invalid[2];cudaMemcpy(invalid,compact_counts_g,2*sizeof(unsigned long long),cudaMemcpyDeviceToHost);
        if(invalid[0]){double details[4];cudaMemcpy(details,compact_max_g,4*sizeof(double),cudaMemcpyDeviceToHost);const int n=int(invalid[1]);const int k=n/gd.ijcells;const int rem=n-k*gd.ijcells;const int j=rem/gd.icells;const int i=rem-j*gd.icells;throw std::runtime_error("AMD_TIMESTEP_INVALID variable="+name+" sgs="+std::to_string(details[1])+" molecular="+std::to_string(details[2])+" metric="+std::to_string(details[3])+" i="+std::to_string(i)+" j="+std::to_string(j)+" k="+std::to_string(k)+" rank="+std::to_string(master.get_mpiid()));}
    };
    check("momentum",fields.sd.at("evisc")->fld_g,fields.visc);
    for(const auto& item:fields.sp)check(item.first,swamd_scalar?fields.sd.at(scalar_coeff.at(item.first))->fld_g.data():nullptr,item.second->visc);
    double maximum;cudaMemcpy(&maximum,compact_max_g,sizeof(double),cudaMemcpyDeviceToHost);master.max(&maximum,1);return maximum*dt;
}

template<typename TF>
unsigned long Diff_amd2<TF>::get_time_limit(const unsigned long idt,const double dt)
{
    const double maximum=get_dn(1.);return idt*dnmax/(std::max(maximum,double(Constants::dsmall))*dt);
}

#ifdef FLOAT_SINGLE
template class Diff_amd2<float>;
#else
template class Diff_amd2<double>;
#endif
