/*
 * MicroHH AMD2 diffusion model -- strict host implementation.
 */
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize ("no-fast-math", "no-finite-math-only", "no-associative-math", "fp-contract=off")
#elif defined(__clang__)
#pragma clang fp reassociate(off)
#pragma clang fp contract(off)
#elif !defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER)
#error "AMD2 requires an explicit strict floating-point policy for this compiler"
#endif

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

#include <netcdf.h>

#include "boundary.h"
#include "constants.h"
#include "diff_amd2.h"
#include "diff_kernels.h"
#include "diff_amd2_math.h"
#include "fields.h"
#include "grid.h"
#include "input.h"
#include "master.h"
#include "stats.h"
#include "thermo.h"

namespace
{
    namespace dk = Diff_kernels;

    inline double sub(const double a, const double b) { return a-b; }

    template<typename TF>
    double centered_x(const TF* f, const int n, const double dxi)
    {
        const double east = sub(double(f[n+1]), double(f[n]));
        const double west = sub(double(f[n]), double(f[n-1]));
        return (east+west)*(0.5*dxi);
    }

    template<typename TF>
    double centered_y(const TF* f, const int n, const int jj, const double dyi)
    {
        const double north = sub(double(f[n+jj]), double(f[n]));
        const double south = sub(double(f[n]), double(f[n-jj]));
        return (north+south)*(0.5*dyi);
    }

    template<typename TF>
    double centered_z(const TF* f, const int n, const int k, const int kk,
                      const double alpha, const TF* dzhi)
    {
        const double lower = sub(double(f[n]), double(f[n-kk]))*double(dzhi[k]);
        const double upper = sub(double(f[n+kk]), double(f[n]))*double(dzhi[k+1]);
        const double lower_weighted = (1.-alpha)*lower;
        const double upper_weighted = alpha*upper;
        return lower_weighted+upper_weighted;
    }

    template<typename TF>
    double one_sided_z(const TF* f, const int n, const int kk,
                       const TF* z, const int k0, const int k1, const int k2)
    {
        const double x0=z[k0], x1=z[k1], x2=z[k2];
        const double w0=(2.*x0-x1-x2)/((x0-x1)*(x0-x2));
        const double w1=(x0-x2)/((x1-x0)*(x1-x2));
        const double w2=(x0-x1)/((x2-x0)*(x2-x1));
        return w0*double(f[n]) + w1*double(f[n+(k1-k0)*kk]) + w2*double(f[n+(k2-k0)*kk]);
    }

    template<typename TF>
    bool representable(const double value)
    {
        return std::isfinite(value) && std::abs(value) <= double(std::numeric_limits<TF>::max());
    }

    template<typename TF, Surface_model surface_model>
    void diff_c_molecular_les(
            TF* at, const TF* a, const TF* dzi, const TF* dzhi,
            const TF dxidxi, const TF dyidyi, const TF* fluxbot, const TF* fluxtop,
            const TF* rhoref, const TF* rhorefh, const TF visc,
            const int istart, const int iend, const int jstart, const int jend,
            const int kstart, const int kend, const int jj, const int kk)
    {
        constexpr int ko = surface_model == Surface_model::Disabled ? 0 : 1;
        if constexpr (surface_model == Surface_model::Enabled)
        {
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ij=i+j*jj, n=ij+kstart*kk;
                    at[n] += visc*((a[n+1]-TF(2)*a[n]+a[n-1])*dxidxi
                                +(a[n+jj]-TF(2)*a[n]+a[n-jj])*dyidyi)
                           +(rhorefh[kstart+1]*visc*(a[n+kk]-a[n])*dzhi[kstart+1]
                            +rhorefh[kstart]*fluxbot[ij])/rhoref[kstart]*dzi[kstart];
                }
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ij=i+j*jj, n=ij+(kend-1)*kk;
                    at[n] += visc*((a[n+1]-TF(2)*a[n]+a[n-1])*dxidxi
                                +(a[n+jj]-TF(2)*a[n]+a[n-jj])*dyidyi)
                           +(-rhorefh[kend]*fluxtop[ij]
                             -rhorefh[kend-1]*visc*(a[n]-a[n-kk])*dzhi[kend-1])
                             /rhoref[kend-1]*dzi[kend-1];
                }
        }
        for (int k=kstart+ko; k<kend-ko; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int n=i+j*jj+k*kk;
                    at[n] += visc*((a[n+1]-TF(2)*a[n]+a[n-1])*dxidxi
                                +(a[n+jj]-TF(2)*a[n]+a[n-jj])*dyidyi)
                           +(rhorefh[k+1]*visc*(a[n+kk]-a[n])*dzhi[k+1]
                            -rhorefh[k]*visc*(a[n]-a[n-kk])*dzhi[k])
                             /rhoref[k]*dzi[k];
                }
    }

    template<typename TF, Surface_model surface_model>
    void flux_c_molecular_les(TF* out, const TF* a, const TF* dzhi, const TF visc,
                              const int istart, const int iend, const int jstart, const int jend,
                              const int kstart, const int kend, const int jj, const int kk)
    {
        constexpr int ko=surface_model == Surface_model::Disabled ? 0 : 1;
        for (int k=kstart+ko; k<=kend-ko; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int n=i+j*jj+k*kk;
                    out[n]=-visc*(a[n]-a[n-kk])*dzhi[k];
                }
    }
}

template<typename TF>
Diff_amd2<TF>::Diff_amd2(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
                         Boundary<TF>& boundaryin, Input& inputin) :
    Diff<TF>(masterin, gridin, fieldsin, boundaryin, inputin),
    boundary_cyclic(masterin, gridin), input(inputin)
{
    dnmax=input.get_item<double>("diff", "dnmax", "", 0.4);
    camd[0]=input.get_item<double>("diff", "camd_x", "", 1./3.);
    camd[1]=input.get_item<double>("diff", "camd_y", "", 1./3.);
    camd[2]=input.get_item<double>("diff", "camd_z", "", 1./3.);
    swamd_buoy=input.get_item<bool>("diff", "swamd_buoy", "", true);
    swamd_scalar=input.get_item<bool>("diff", "swamd_scalar", "", true);
    amd_max=input.get_item<double>("diff", "amd_max", "", 0.);
    if (grid.get_spatial_order()!=Grid_order::Second)
        throw std::runtime_error("swdiff=amd2 requires swspatialorder=2");
}

template<typename TF>
std::string Diff_amd2<TF>::base32_token(const std::string& value)
{
    return AMD2_math::base32_token(value);
}

template<typename TF>
void Diff_amd2<TF>::register_fields()
{
    auto& gd=grid.get_grid_data();
    std::vector<std::string> proposed={"evisc", "nu_s", "nu_b"};
    for (const auto& item : fields.sp)
        if (swamd_scalar) proposed.push_back("ediff_"+item.first);

    const auto namespaces=[&](const std::string& name)
    {
        std::vector<std::string> found;
        if (fields.a.count(name)) found.push_back("a");
        if (fields.ap.count(name)) found.push_back("ap");
        if (fields.at.count(name)) found.push_back("at");
        if (fields.mp.count(name)) found.push_back("mp");
        if (fields.mt.count(name)) found.push_back("mt");
        if (fields.sp.count(name)) found.push_back("sp");
        if (fields.st.count(name)) found.push_back("st");
        if (fields.sd.count(name)) found.push_back("sd");
        return found;
    };
    std::set<std::string> unique;
    for (const auto& name : proposed)
    {
        const auto found=namespaces(name);
        if (!unique.insert(name).second || !found.empty())
        {
            std::ostringstream msg;
            msg << "AMD field collision for proposed field '" << name << "'";
            for (const auto& ns : found) msg << " " << ns;
            throw std::runtime_error(msg.str());
        }
    }

    std::set<std::string> tokens,netcdf_names;
    static const char* scalar_statuses[]={"zero_gradient","invalid_gradient","invalid_contraction","clipped","storage_overflow","positive","capped"};
    for (const auto& item : fields.sp)
    {
        const std::string token=base32_token(item.first);
        const std::string profile="amd_ediff_"+token;
        if (!tokens.insert(token).second)
            throw std::runtime_error("AMD Base32 scalar token collision: "+token);
        std::vector<std::string> encoded={profile,profile+"_2","amd_scalar_name_"+token,"amd_"+token+"_max"};
        for(const char* status:scalar_statuses){encoded.push_back("amd_"+token+"_"+status+"_count");encoded.push_back("amd_"+token+"_"+status+"_fraction");}
        for(const auto& name:encoded)
        {
            if(name.size()>NC_MAX_NAME)throw std::runtime_error("AMD encoded name exceeds NC_MAX_NAME for scalar '"+item.first+"': "+name);
            if(!netcdf_names.insert(name).second)throw std::runtime_error("AMD encoded NetCDF name collision: "+name);
        }
        scalar_token[item.first]=token;
        if (swamd_scalar) scalar_coeff[item.first]="ediff_"+item.first;
    }

    fields.init_diagnostic_field("evisc", "AMD eddy viscosity", "m2 s-1", "default", gd.sloc);
    fields.init_diagnostic_field("nu_s", "AMD shear contribution", "m2 s-1", "default", gd.sloc);
    fields.init_diagnostic_field("nu_b", "AMD buoyancy contribution", "m2 s-1", "default", gd.sloc);
    for (const auto& item : scalar_coeff)
        fields.init_diagnostic_field(item.second, "AMD diffusivity for "+item.first, "m2 s-1", "default", gd.sloc);
}

template<typename TF>
void Diff_amd2<TF>::validate_configuration(const Thermo<TF>& thermoin)
{
    thermo=&thermoin;
    if (!std::isfinite(dnmax) || dnmax<=0.) throw std::runtime_error("AMD dnmax must be finite and positive");
    for (int n=0; n<3; ++n)
        if (!std::isfinite(camd[n]) || camd[n]<=0.) throw std::runtime_error("AMD camd_x/y/z must be finite and positive");
    if (!std::isfinite(amd_max)) throw std::runtime_error("AMD amd_max must be finite");
    if (amd_max>0. && !representable<TF>(amd_max)) throw std::runtime_error("AMD amd_max is not representable in model precision");
    if (!std::isfinite(double(fields.visc)) || fields.visc<TF(0)) throw std::runtime_error("AMD molecular momentum viscosity is invalid");
    for (const auto& item : fields.sp)
        if (!std::isfinite(double(item.second->visc)) || item.second->visc<TF(0))
            throw std::runtime_error("AMD molecular diffusivity is invalid for scalar '"+item.first+"'");
    const auto& gd=grid.get_grid_data();
    const bool ustar=boundary.get_momentum_bcbot()==Boundary_type::Ustar_type
                  || boundary.get_momentum_bctop()==Boundary_type::Ustar_type;
    if (ustar && gd.kmax<3) throw std::runtime_error("AMD Ustar boundary requires at least three physical vertical levels");
    if (swamd_buoy && thermo->get_switch()!=Thermo_type::Disabled
            && (!thermo->check_field_exists("b") || !thermo->check_field_exists("b_h")))
        throw std::runtime_error("AMD buoyancy requires thermodynamic fields b and b_h");

    const auto geometry=thermo->get_buoyancy_geometry();
    std::ostringstream msg;
    msg << "AMD2: C=(" << camd[0] << "," << camd[1] << "," << camd[2]
        << "), buoyancy=" << swamd_buoy << ", scalars=" << swamd_scalar
        << ", cap=" << (amd_max>0. ? std::to_string(amd_max) : std::string("disabled"))
        << ", strict-fp=on, d=(" << geometry.force_direction[0] << ","
        << geometry.force_direction[1] << "," << geometry.force_direction[2]
        << "); no fourth-order/IB/wall damping";
    master.print_message(msg.str());
}

template<typename TF> void Diff_amd2<TF>::init() { boundary_cyclic.init(); }

template<typename TF>
void Diff_amd2<TF>::fill_coefficient_halos(Field3d<TF>& fld)
{
    auto& gd=grid.get_grid_data();
    boundary_cyclic.exec(fld.fld.data());
    for (int j=0; j<gd.jcells; ++j)
        for (int i=0; i<gd.icells; ++i)
        {
            const int b=i+j*gd.icells+gd.kstart*gd.ijcells;
            const int t=i+j*gd.icells+(gd.kend-1)*gd.ijcells;
            fld.fld[b-gd.ijcells]=fld.fld[b];
            fld.fld[t+gd.ijcells]=fld.fld[t];
        }
}

template<typename TF>
void Diff_amd2<TF>::calculate_coefficients(Thermo<TF>& thermoin, const Mode mode)
{
    auto& gd=grid.get_grid_data();
    auto& evisc=*fields.sd.at("evisc");
    auto& nus=*fields.sd.at("nu_s");
    auto& nub=*fields.sd.at("nu_b");
    const TF* u=fields.mp.at("u")->fld.data();
    const TF* v=fields.mp.at("v")->fld.data();
    const TF* w=fields.mp.at("w")->fld.data();
    std::shared_ptr<Field3d<TF>> b, bh;
    const bool use_buoy=swamd_buoy && thermoin.get_switch()!=Thermo_type::Disabled;
    if (use_buoy)
    {
        b=fields.get_tmp(); bh=fields.get_tmp();
        thermoin.get_thermo_field(*b, "b", true, false);
        thermoin.get_thermo_field(*bh, "b_h", false, false);
    }
    const auto geometry=thermoin.get_buoyancy_geometry();
    const double d[3]={double(geometry.force_direction[0]),double(geometry.force_direction[1]),double(geometry.force_direction[2])};
    const double G[3]={double(geometry.background_gradient[0]),double(geometry.background_gradient[1]),double(geometry.background_gradient[2])};
    const double Mxy[2]={camd[0]*double(gd.dx)*double(gd.dx),camd[1]*double(gd.dy)*double(gd.dy)};
    const bool ustar_bot=boundary.get_momentum_bcbot()==Boundary_type::Ustar_type;
    const bool ustar_top=boundary.get_momentum_bctop()==Boundary_type::Ustar_type;
    Aggregates pass;
    for (const auto& item : scalar_coeff) pass.scalar[item.first]=std::vector<std::uint64_t>(Scalar_status_count,0);

    for (int k=gd.kstart; k<gd.kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int n=i+j*gd.icells+k*gd.ijcells;
                evisc.fld[n]=TF(0); nus.fld[n]=TF(0); nub.fld[n]=TF(0);
                for (const auto& item : scalar_coeff) fields.sd.at(item.second)->fld[n]=TF(0);
                Mom_status status=Mom_invalid_velocity;
                bool scalar_prerequisite=false;
                const double alpha=(double(gd.z[k])-double(gd.zh[k]))/(double(gd.zh[k+1])-double(gd.zh[k]));
                double A[3][3];
                A[0][0]=sub(double(u[n+1]),double(u[n]))*double(gd.dxi);
                A[1][1]=sub(double(v[n+gd.icells]),double(v[n]))*double(gd.dyi);
                A[2][2]=sub(double(w[n+gd.ijcells]),double(w[n]))*double(gd.dzi[k]);
                const double udy0=sub(double(u[n+gd.icells]),double(u[n-gd.icells]));
                const double udy1=sub(double(u[n+1+gd.icells]),double(u[n+1-gd.icells]));
                A[0][1]=(udy0+udy1)*(0.25*double(gd.dyi));
                const double vdx0=sub(double(v[n+1]),double(v[n-1]));
                const double vdx1=sub(double(v[n+gd.icells+1]),double(v[n+gd.icells-1]));
                A[1][0]=(vdx0+vdx1)*(0.25*double(gd.dxi));
                auto uz=[&](const TF* f, const int q)
                {
                    if (k==gd.kstart && ustar_bot) return one_sided_z(f,q,gd.ijcells,gd.z.data(),k,k+1,k+2);
                    if (k==gd.kend-1 && ustar_top) return one_sided_z(f,q,gd.ijcells,gd.z.data(),k,k-1,k-2);
                    return centered_z(f,q,k,gd.ijcells,alpha,gd.dzhi.data());
                };
                const double uz0=uz(u,n), uz1=uz(u,n+1); A[0][2]=0.5*(uz0+uz1);
                const double vz0=uz(v,n), vz1=uz(v,n+gd.icells); A[1][2]=0.5*(vz0+vz1);
                const double wx0=centered_x(w,n,double(gd.dxi));
                const double wx1=centered_x(w,n+gd.ijcells,double(gd.dxi));
                A[2][0]=(1.-alpha)*wx0+alpha*wx1;
                const double wy0=centered_y(w,n,gd.icells,double(gd.dyi));
                const double wy1=centered_y(w,n+gd.ijcells,gd.icells,double(gd.dyi));
                A[2][1]=(1.-alpha)*wy0+alpha*wy1;

                double a=0.; bool finite_A=std::isfinite(alpha);
                for (int ii=0; ii<3; ++ii) for (int jj=0; jj<3; ++jj)
                { finite_A=finite_A&&std::isfinite(A[ii][jj]); a=std::max(a,std::abs(A[ii][jj])); }
                if (!finite_A || !std::isfinite(a)) status=Mom_invalid_velocity;
                else if (a==0.) status=Mom_zero_velocity;
                else
                {
                    double Ah[3][3], Sh[3][3], Du=0.;
                    for (int ii=0; ii<3; ++ii) for (int jj=0; jj<3; ++jj)
                    { Ah[ii][jj]=A[ii][jj]/a; Du+=Ah[ii][jj]*Ah[ii][jj]; }
                    for (int ii=0; ii<3; ++ii) for (int jj=0; jj<3; ++jj) Sh[ii][jj]=0.5*(Ah[ii][jj]+Ah[jj][ii]);
                    if (!std::isfinite(Du) || !(Du>0.)) status=Mom_invalid_velocity;
                    else
                    {
                        scalar_prerequisite=true;
                        const double M[3]={Mxy[0],Mxy[1],camd[2]*double(gd.dz[k])*double(gd.dz[k])};
                        const double shear=AMD2_math::shear_contraction(Ah,Sh,M);
                        const double nu_s=-a*shear/Du;
                        double nu_b=0.; bool zero_b=false;
                        if (!std::isfinite(shear) || !std::isfinite(nu_s)) status=Mom_invalid_shear;
                        else
                        {
                            bool buoy_ok=true;
                            if (use_buoy)
                            {
                                double q[3];
                                q[0]=centered_x(b->fld.data(),n,double(gd.dxi))+G[0];
                                q[1]=centered_y(b->fld.data(),n,gd.icells,double(gd.dyi))+G[1];
                                q[2]=sub(double(bh->fld[n+gd.ijcells]),double(bh->fld[n]))*double(gd.dzi[k])+G[2];
                                double c=0.; for (int kk=0; kk<3; ++kk) { buoy_ok=buoy_ok&&std::isfinite(q[kk]); c=std::max(c,std::abs(q[kk])); }
                                if (!std::isfinite(c)) buoy_ok=false;
                                else if (c==0.) zero_b=true;
                                else
                                {
                                    double qhat[3]; for (int kk=0;kk<3;++kk) qhat[kk]=q[kk]/c;
                                    const double contraction=AMD2_math::buoyancy_contraction(Ah,qhat,d,M);
                                    nu_b=(c/a)*contraction/Du;
                                    buoy_ok=buoy_ok&&std::isfinite(contraction)&&std::isfinite(nu_b);
                                }
                            }
                            if (!buoy_ok) status=Mom_invalid_buoyancy;
                            else if (!representable<TF>(nu_s) || !representable<TF>(nu_b)) status=Mom_diagnostic_storage_overflow;
                            else
                            {
                                const double raw=nu_s+nu_b;
                                if (!std::isfinite(raw)) status=Mom_invalid_sum;
                                else if (raw<=0.) { nus.fld[n]=TF(nu_s); nub.fld[n]=TF(nu_b); status=Mom_clipped; }
                                else
                                {
                                    const double stored=(amd_max>0. && raw>amd_max) ? amd_max : raw;
                                    if (!representable<TF>(stored)) status=Mom_storage_overflow;
                                    else
                                    {
                                        nus.fld[n]=TF(nu_s); nub.fld[n]=TF(nu_b); evisc.fld[n]=TF(stored);
                                        status=(stored!=raw) ? Mom_capped : Mom_positive;
                                        pass.max_evisc=std::max(pass.max_evisc,stored);
                                    }
                                }
                            }
                            if (zero_b && mode==Mode::Operational) ++pass.zero_buoyancy_gradient;
                        }

                        if (scalar_prerequisite)
                            for (const auto& item : scalar_coeff)
                            {
                                auto& counts=pass.scalar[item.first];
                                auto& ediff=*fields.sd.at(item.second);
                                const TF* phi=fields.sp.at(item.first)->fld.data();
                                double grad[3]={centered_x(phi,n,double(gd.dxi)),centered_y(phi,n,gd.icells,double(gd.dyi)),centered_z(phi,n,k,gd.ijcells,alpha,gd.dzhi.data())};
                                double g=0.; bool ok=true; for (double x:grad) { ok=ok&&std::isfinite(x); g=std::max(g,std::abs(x)); }
                                Scalar_status ss=Scalar_invalid_gradient;
                                if (!ok || !std::isfinite(g)) ss=Scalar_invalid_gradient;
                                else if (g==0.) ss=Scalar_zero_gradient;
                                else
                                {
                                    double p[3], D=0.; for (int ii=0;ii<3;++ii) { p[ii]=grad[ii]/g; D+=p[ii]*p[ii]; }
                                    if (!std::isfinite(D) || !(D>0.)) ss=Scalar_invalid_gradient;
                                    else
                                    {
                                        const double sum=AMD2_math::scalar_contraction(Ah,p,M);
                                        const double raw=-a*sum/D;
                                        if (!std::isfinite(sum) || !std::isfinite(raw)) ss=Scalar_invalid_contraction;
                                        else if (raw<=0.) ss=Scalar_clipped;
                                        else
                                        {
                                            const double stored=(amd_max>0.&&raw>amd_max)?amd_max:raw;
                                            if (!representable<TF>(stored)) ss=Scalar_storage_overflow;
                                            else { ediff.fld[n]=TF(stored); ss=stored!=raw?Scalar_capped:Scalar_positive; pass.max_scalar[item.first]=std::max(pass.max_scalar[item.first],stored); }
                                        }
                                    }
                                }
                                if (mode==Mode::Operational) ++counts[ss];
                            }
                    }
                }
                if (mode==Mode::Operational) { ++pass.cells_evaluated; ++pass.mom[status]; }
            }

    fill_coefficient_halos(evisc); fill_coefficient_halos(nus); fill_coefficient_halos(nub);
    for (const auto& item : scalar_coeff) fill_coefficient_halos(*fields.sd.at(item.second));
    if (use_buoy) { fields.release_tmp(b); fields.release_tmp(bh); }
    if (mode==Mode::Operational) { add_aggregates(active,pass); add_aggregates(cumulative,pass); }
}

template<typename TF>
void Diff_amd2<TF>::add_aggregates(Aggregates& dst, const Aggregates& src)
{
    const auto checked_add=[](std::uint64_t& destination,const std::uint64_t value)
    {
        if (value>std::numeric_limits<std::uint64_t>::max()-destination)
            throw std::overflow_error("AMD uint64 diagnostic counter overflow");
        destination+=value;
    };
    checked_add(dst.cells_evaluated,src.cells_evaluated);
    for (int n=0;n<Mom_status_count;++n) checked_add(dst.mom[n],src.mom[n]);
    checked_add(dst.zero_buoyancy_gradient,src.zero_buoyancy_gradient);
    dst.max_evisc=std::max(dst.max_evisc,src.max_evisc);
    for (const auto& item:src.scalar)
    { auto& v=dst.scalar[item.first]; if (v.empty()) v.resize(Scalar_status_count); for (int n=0;n<Scalar_status_count;++n) checked_add(v[n],item.second[n]); }
    for (const auto& item:src.max_scalar) dst.max_scalar[item.first]=std::max(dst.max_scalar[item.first],item.second);
}

template<typename TF>
void Diff_amd2<TF>::create(Stats<TF>& stats, const bool cold_start)
{
    if (cold_start || !stats.get_switch()) return;
    const std::string group="default";
    stats.add_profs(*fields.sd.at("evisc"),"z",{"mean","2"},group);
    stats.add_profs(*fields.sd.at("nu_s"),"z",{"mean","2"},group);
    stats.add_profs(*fields.sd.at("nu_b"),"z",{"mean","2"},group);
    for (const auto& item:scalar_coeff)
    {
        const std::string alias="amd_ediff_"+scalar_token.at(item.first);
        stats.add_prof(alias,"AMD diffusivity for "+item.first,"m2 s-1","z",group);
        stats.add_prof(alias+"_2","AMD diffusivity squared for "+item.first,"m4 s-2","z",group);
    }
    stats.add_tendency(*fields.mt.at("u"),"z","diff","Diffusion");
    stats.add_tendency(*fields.mt.at("v"),"z","diff","Diffusion");
    stats.add_tendency(*fields.mt.at("w"),"zh","diff","Diffusion");
    for (const auto& item:fields.st) stats.add_tendency(*item.second,"z","diff","Diffusion");

    static const char* mom_names[Mom_status_count]={
        "invalid_velocity","zero_velocity","invalid_shear","invalid_buoyancy",
        "diagnostic_storage_overflow","invalid_sum","clipped","storage_overflow","positive","capped"};
    static const char* scalar_names[Scalar_status_count]={
        "zero_gradient","invalid_gradient","invalid_contraction","clipped","storage_overflow","positive","capped"};
    stats.add_time_series_u64("amd_cells_evaluated","Domain-wide AMD operational cells evaluated","1",group);
    for (int n=0;n<Mom_status_count;++n)
    {
        const std::string base="amd_mom_"+std::string(mom_names[n]);
        stats.add_time_series_u64(base+"_count","Domain-wide AMD momentum "+std::string(mom_names[n])+" count","1",group);
        stats.add_time_series(base+"_fraction","AMD momentum "+std::string(mom_names[n])+" fraction","1",group);
    }
    stats.add_time_series_u64("amd_mom_zero_buoyancy_gradient_count","Domain-wide zero buoyancy-gradient count","1",group);
    stats.add_time_series("amd_mom_zero_buoyancy_gradient_fraction","Zero buoyancy-gradient fraction","1",group);
    stats.add_time_series("amd_mom_max","Maximum AMD momentum viscosity","m2 s-1",group);
    for (const auto& item:scalar_coeff)
    {
        const std::string& token=scalar_token.at(item.first);
        const std::string base="amd_"+token;
        for (int n=0;n<Scalar_status_count;++n)
        {
            const std::string status=base+"_"+scalar_names[n];
            stats.add_time_series_u64(status+"_count","Domain-wide AMD scalar "+std::string(scalar_names[n])+" count","1",group);
            stats.add_time_series(status+"_fraction","AMD scalar "+std::string(scalar_names[n])+" fraction","1",group);
        }
        stats.add_time_series(base+"_max","Maximum AMD diffusivity for "+item.first,"m2 s-1",group);
        stats.add_global_attribute("amd_scalar_name_"+token,item.first);
    }
    stats.add_global_attribute("amd_model","amd2");
    stats.add_global_attribute("amd_formulation_version","AMD2 normalized-generalized-buoyancy-v1");
    stats.add_global_attribute("amd_precision_policy","strict; reassociation, contraction, finite-only assumptions and fast approximations disabled");
    stats.add_global_attribute("amd_storage_precision",sizeof(TF)==sizeof(float)?"float32":"float64");
    stats.add_global_attribute("amd_contraction_precision","float64");
    stats.add_global_attribute("amd_fma_enabled","false");
    stats.add_global_attribute("amd_counts_mask_conditioned","false");
    stats.add_global_attribute("amd_cumulative_counter_scope","run_segment_not_restart_persistent");
    stats.add_global_attribute("amd_counter_scope","domain-wide operational RK-stage cells; repeated identically in every mask file");
    stats.add_global_attribute("amd_surface_velocity_gradient","resolved_three_point_for_ustar");
    stats.add_global_attribute("amd_surface_gradient_policy","resolved ghosts; three-point stretched-grid one-sided derivative for momentum Ustar boundaries");
    stats.add_global_attribute("amd_unsupported_features","fourth-order grids, active immersed boundaries, wall damping");
    stats.add_global_attribute("amd_budget_semantics","Budget_2 historical Smagorinsky LES helper; *_diff omits parts of the molecular/anelastic operator and is not an exact total-diffusion closure");
    stats.add_global_attribute("amd_spatial_order",2.);
    stats.add_global_attribute("amd_immersed_boundary_supported","false");
    stats.add_global_attribute("amd_cx",camd[0]); stats.add_global_attribute("amd_cy",camd[1]); stats.add_global_attribute("amd_cz",camd[2]);
    stats.add_global_attribute("amd_max",amd_max); stats.add_global_attribute("amd_cap_enabled",amd_max>0.?1.:0.);
    stats.add_global_attribute("amd_buoyancy_enabled",swamd_buoy?1.:0.); stats.add_global_attribute("amd_scalar_enabled",swamd_scalar?1.:0.);
    if (thermo)
    {
        const auto g=thermo->get_buoyancy_geometry();
        for (int n=0;n<3;++n) { stats.add_global_attribute("amd_buoyancy_direction_"+std::to_string(n+1),double(g.force_direction[n])); stats.add_global_attribute("amd_background_gradient_"+std::to_string(n+1),double(g.background_gradient[n])); }
    }
}

#ifndef USECUDA
template<typename TF> void Diff_amd2<TF>::exec_viscosity(Stats<TF>&, Thermo<TF>& th) { calculate_coefficients(th,Mode::Operational); }

template<typename TF>
void Diff_amd2<TF>::exec(Stats<TF>& stats)
{
    auto& gd=grid.get_grid_data();
    auto wrapper=[&]<Surface_model sm>()
    {
        const TF* ev=fields.sd.at("evisc")->fld.data();
        dk::diff_u<TF,sm>(fields.mt.at("u")->fld.data(),fields.mp.at("u")->fld.data(),fields.mp.at("v")->fld.data(),fields.mp.at("w")->fld.data(),gd.dzi.data(),gd.dzhi.data(),1./gd.dx,1./gd.dy,ev,fields.mp.at("u")->flux_bot.data(),fields.mp.at("u")->flux_top.data(),fields.rhoref.data(),fields.rhorefh.data(),fields.visc,gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
        dk::diff_v<TF,sm>(fields.mt.at("v")->fld.data(),fields.mp.at("u")->fld.data(),fields.mp.at("v")->fld.data(),fields.mp.at("w")->fld.data(),gd.dzi.data(),gd.dzhi.data(),1./gd.dx,1./gd.dy,ev,fields.mp.at("v")->flux_bot.data(),fields.mp.at("v")->flux_top.data(),fields.rhoref.data(),fields.rhorefh.data(),fields.visc,gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
        dk::diff_w<TF>(fields.mt.at("w")->fld.data(),fields.mp.at("u")->fld.data(),fields.mp.at("v")->fld.data(),fields.mp.at("w")->fld.data(),gd.dzi.data(),gd.dzhi.data(),1./gd.dx,1./gd.dy,ev,fields.rhoref.data(),fields.rhorefh.data(),fields.visc,gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
        for (const auto& item:fields.st)
        {
            auto& p=*fields.sp.at(item.first);
            if (swamd_scalar)
                dk::diff_c<TF,sm>(item.second->fld.data(),p.fld.data(),gd.dzi.data(),gd.dzhi.data(),1./(gd.dx*gd.dx),1./(gd.dy*gd.dy),fields.sd.at(scalar_coeff.at(item.first))->fld.data(),p.flux_bot.data(),p.flux_top.data(),fields.rhoref.data(),fields.rhorefh.data(),TF(1),p.visc,gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
            else
                diff_c_molecular_les<TF,sm>(item.second->fld.data(),p.fld.data(),gd.dzi.data(),gd.dzhi.data(),1./(gd.dx*gd.dx),1./(gd.dy*gd.dy),p.flux_bot.data(),p.flux_top.data(),fields.rhoref.data(),fields.rhorefh.data(),p.visc,gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
        }
    };
    if (boundary.get_switch()=="default") wrapper.template operator()<Surface_model::Disabled>(); else wrapper.template operator()<Surface_model::Enabled>();
    for (const auto& n:{"u","v","w"}) stats.calc_tend(*fields.mt.at(n),"diff");
    for (const auto& item:fields.st) stats.calc_tend(*item.second,"diff");
}

#endif

template<typename TF>
void Diff_amd2<TF>::diff_flux(Field3d<TF>& out, const Field3d<TF>& in)
{
    auto& gd=grid.get_grid_data();
    auto wrapper=[&]<Surface_model sm>()
    {
        if (in.loc[0]==1) dk::calc_diff_flux_u<TF,sm>(out.fld.data(),in.fld.data(),fields.mp.at("w")->fld.data(),fields.sd.at("evisc")->fld.data(),gd.dxi,gd.dzhi.data(),fields.visc,gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
        else if (in.loc[1]==1) dk::calc_diff_flux_v<TF,sm>(out.fld.data(),in.fld.data(),fields.mp.at("w")->fld.data(),fields.sd.at("evisc")->fld.data(),gd.dyi,gd.dzhi.data(),fields.visc,gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
        else if (fields.mp.count(in.name)) dk::calc_diff_flux_c<TF,sm>(out.fld.data(),in.fld.data(),fields.sd.at("evisc")->fld.data(),gd.dzhi.data(),TF(1),fields.visc,gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
        else if (swamd_scalar && scalar_coeff.count(in.name)) dk::calc_diff_flux_c<TF,sm>(out.fld.data(),in.fld.data(),fields.sd.at(scalar_coeff.at(in.name))->fld.data(),gd.dzhi.data(),TF(1),in.visc,gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
        else flux_c_molecular_les<TF,sm>(out.fld.data(),in.fld.data(),gd.dzhi.data(),in.visc,gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.kend,gd.icells,gd.ijcells);
    };
    const bool surface=boundary.get_switch()!="default";
    if (surface)
    {
        dk::calc_diff_flux_bc(out.fld.data(),in.flux_bot.data(),gd.istart,gd.iend,gd.jstart,gd.jend,gd.kstart,gd.icells,gd.ijcells);
        dk::calc_diff_flux_bc(out.fld.data(),in.flux_top.data(),gd.istart,gd.iend,gd.jstart,gd.jend,gd.kend,gd.icells,gd.ijcells);
        wrapper.template operator()<Surface_model::Enabled>();
    }
    else wrapper.template operator()<Surface_model::Disabled>();
}

#ifndef USECUDA
template<typename TF>
double Diff_amd2<TF>::timestep_multiplier() const
{
    auto& gd=grid.get_grid_data(); double maximum=0.;
    for (int k=gd.kstart;k<gd.kend;++k) for (int j=gd.jstart;j<gd.jend;++j) for (int i=gd.istart;i<gd.iend;++i)
    {
        const int n=i+j*gd.icells+k*gd.ijcells;
        const double H=1./(double(gd.dx)*double(gd.dx))+1./(double(gd.dy)*double(gd.dy))+double(gd.dzi[k])*double(gd.dzi[k]);
        auto check=[&](const std::string& name,const double sgs,const double molecular)
        { const double value=(sgs+molecular)*H; if (!std::isfinite(sgs)||sgs<0.||!std::isfinite(value)||value<0.) { std::ostringstream m; m<<"AMD_TIMESTEP_INVALID variable="<<name<<" sgs="<<sgs<<" molecular="<<molecular<<" metric="<<H<<" i="<<i<<" j="<<j<<" k="<<k<<" rank="<<master.get_mpiid(); throw std::runtime_error(m.str()); } maximum=std::max(maximum,value); };
        check("momentum",double(fields.sd.at("evisc")->fld[n]),double(fields.visc));
        for (const auto& item:fields.sp) check(item.first,swamd_scalar?double(fields.sd.at(scalar_coeff.at(item.first))->fld[n]):0.,double(item.second->visc));
    }
    master.max(&maximum,1); return maximum;
}
template<typename TF> unsigned long Diff_amd2<TF>::get_time_limit(const unsigned long idt,const double dt) { return idt*dnmax/(std::max(timestep_multiplier(),double(Constants::dsmall))*dt); }
template<typename TF> double Diff_amd2<TF>::get_dn(const double dt) { return timestep_multiplier()*dt; }
#endif

template<typename TF>
void Diff_amd2<TF>::prepare_stats()
{
    if (pending_valid) throw std::runtime_error("AMD statistics snapshot was not consumed before the next output event");
    pending=std::move(active); active=Aggregates{}; pending_valid=true;
}

template<typename TF>
void Diff_amd2<TF>::exec_stats(Stats<TF>& stats,Thermo<TF>& th)
{
    calculate_coefficients(th,Mode::Diagnostic);
    stats.calc_stats("evisc",*fields.sd.at("evisc"),TF(0),TF(0));
    stats.calc_stats("nu_s",*fields.sd.at("nu_s"),TF(0),TF(0));
    stats.calc_stats("nu_b",*fields.sd.at("nu_b"),TF(0),TF(0));
    for (const auto& item:scalar_coeff)
    {
        const std::string alias="amd_ediff_"+scalar_token.at(item.first);
        stats.calc_stats(alias,*fields.sd.at(item.second),TF(0),TF(0));
    }
    if (pending_valid)
    {
        static const char* mom_names[Mom_status_count]={
            "invalid_velocity","zero_velocity","invalid_shear","invalid_buoyancy",
            "diagnostic_storage_overflow","invalid_sum","clipped","storage_overflow","positive","capped"};
        static const char* scalar_names[Scalar_status_count]={
            "zero_gradient","invalid_gradient","invalid_contraction","clipped","storage_overflow","positive","capped"};
        master.sum(&pending.cells_evaluated,1);
        stats.set_time_series_u64("amd_cells_evaluated",pending.cells_evaluated);
        master.sum(pending.mom,Mom_status_count);
        master.sum(&pending.zero_buoyancy_gradient,1);
        for (int n=0;n<Mom_status_count;++n)
        {
            const std::string base="amd_mom_"+std::string(mom_names[n]);
            stats.set_time_series_u64(base+"_count",pending.mom[n]);
            stats.set_time_series(base+"_fraction",pending.cells_evaluated?TF(double(pending.mom[n])/double(pending.cells_evaluated)):TF(0));
        }
        stats.set_time_series_u64("amd_mom_zero_buoyancy_gradient_count",pending.zero_buoyancy_gradient);
        stats.set_time_series("amd_mom_zero_buoyancy_gradient_fraction",pending.cells_evaluated?TF(double(pending.zero_buoyancy_gradient)/double(pending.cells_evaluated)):TF(0));
        master.max(&pending.max_evisc,1); stats.set_time_series("amd_mom_max",TF(pending.max_evisc));
        for (auto& item:pending.scalar)
        {
            master.sum(item.second.data(),Scalar_status_count);
            const std::string base="amd_"+scalar_token.at(item.first);
            for (int n=0;n<Scalar_status_count;++n)
            {
                stats.set_time_series_u64(base+"_"+scalar_names[n]+"_count",item.second[n]);
                stats.set_time_series(base+"_"+scalar_names[n]+"_fraction",pending.cells_evaluated?TF(double(item.second[n])/double(pending.cells_evaluated)):TF(0));
            }
            double maximum=pending.max_scalar[item.first]; master.max(&maximum,1); stats.set_time_series(base+"_max",TF(maximum));
        }
    }
    pending_valid=false;
}

template<typename TF>
void Diff_amd2<TF>::finalize_diagnostics()
{
    std::uint64_t evaluated=cumulative.cells_evaluated; master.sum(&evaluated,1);
    master.print_message("AMD2 run-segment cells evaluated: "+std::to_string(evaluated));
    std::uint64_t counts[Mom_status_count]; for (int n=0;n<Mom_status_count;++n) counts[n]=cumulative.mom[n];
    master.sum(counts,Mom_status_count);
    std::ostringstream msg; msg<<"AMD2 run-segment momentum counts:"; for (int n=0;n<Mom_status_count;++n) msg<<" "<<counts[n]; master.print_message(msg.str());
    std::uint64_t zero=cumulative.zero_buoyancy_gradient;master.sum(&zero,1);
    master.print_message("AMD2 run-segment zero-buoyancy-gradient count: "+std::to_string(zero));
    for(const auto& item:cumulative.scalar)
    {
        std::vector<std::uint64_t> scalar_counts=item.second;master.sum(scalar_counts.data(),Scalar_status_count);
        std::ostringstream scalar_msg;scalar_msg<<"AMD2 run-segment scalar counts "<<item.first<<":";for(const auto value:scalar_counts)scalar_msg<<" "<<value;master.print_message(scalar_msg.str());
    }
}

#ifdef FLOAT_SINGLE
template class Diff_amd2<float>;
#else
template class Diff_amd2<double>;
#endif
