#ifndef MICROHHC_DIFF_AMD2_KERNELS_CUH
#define MICROHHC_DIFF_AMD2_KERNELS_CUH

#include "cuda_tiling.h"

namespace Diff_amd2_kernels
{
    template<typename TF, bool surface_model_enabled>
    struct diff_c_molecular_g
    {
        DEFINE_GRID_KERNEL("amd2::diff_c_molecular", surface_model_enabled ? 1 : 0,
                           dim3{32, 2, 2}, true)

        template<typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout gd, int i, int j, int k, Level level,
                TF* __restrict__ at, const TF* __restrict__ a,
                const TF* __restrict__ fluxbot, const TF* __restrict__ fluxtop,
                const TF* __restrict__ dzi, const TF* __restrict__ dzhi,
                const TF* __restrict__ rhorefi, const TF* __restrict__ rhorefh,
                const TF dxidxi, const TF dyidyi, const TF visc)
        {
            const int ii=gd.istride;
            const int jj=gd.jstride;
            const int kk=gd.kstride;
            const int ij=i*ii+j*jj;
            const int ijk=ij+k*kk;
            const TF horizontal=visc*((a[ijk+ii]-TF(2)*a[ijk]+a[ijk-ii])*dxidxi
                                     +(a[ijk+jj]-TF(2)*a[ijk]+a[ijk-jj])*dyidyi);

            if (level.distance_to_start()==0 && surface_model_enabled)
                at[ijk]+=horizontal
                    +(rhorefh[k+1]*visc*(a[ijk+kk]-a[ijk])*dzhi[k+1]
                      +rhorefh[k]*fluxbot[ij])*rhorefi[k]*dzi[k];
            else if (level.distance_to_end()==0 && surface_model_enabled)
                at[ijk]+=horizontal
                    +(-rhorefh[k+1]*fluxtop[ij]
                      -rhorefh[k]*visc*(a[ijk]-a[ijk-kk])*dzhi[k])*rhorefi[k]*dzi[k];
            else
                at[ijk]+=horizontal
                    +(rhorefh[k+1]*visc*(a[ijk+kk]-a[ijk])*dzhi[k+1]
                      -rhorefh[k]*visc*(a[ijk]-a[ijk-kk])*dzhi[k])*rhorefi[k]*dzi[k];
        }
    };
}

#endif
