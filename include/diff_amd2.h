/*
 * MicroHH AMD2 diffusion model.
 *
 * The coefficient calculation intentionally lives in a separate translation
 * unit with a strict floating-point policy; see diff_amd2.cxx.
 */
#ifndef DIFF_AMD2_H
#define DIFF_AMD2_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "boundary_cyclic.h"
#include "diff.h"

template<typename> class Stats;

template<typename TF>
class Diff_amd2 : public Diff<TF>
{
    public:
        Diff_amd2(Master&, Grid<TF>&, Fields<TF>&, Boundary<TF>&, Input&);
        ~Diff_amd2() = default;

        Diffusion_type get_switch() const override { return Diffusion_type::Diff_amd2; }
        void register_fields() override;
        void validate_configuration(const Thermo<TF>&) override;
        void init() override;
        void create(Stats<TF>&, const bool) override;
        void exec(Stats<TF>&) override;
        void exec_viscosity(Stats<TF>&, Thermo<TF>&) override;
        void diff_flux(Field3d<TF>&, const Field3d<TF>&) override;
        void exec_stats(Stats<TF>&, Thermo<TF>&) override;
        unsigned long get_time_limit(unsigned long, double) override;
        double get_dn(double) override;
        void prepare_stats() override;
        void finalize_diagnostics() override;

        #ifdef USECUDA
        void prepare_device(Boundary<TF>&) override;
        void clear_device() override;
        cuda_vector<unsigned long long> compact_counts_g;
        cuda_vector<double> compact_max_g;
        #endif

    private:
        using Diff<TF>::master;
        using Diff<TF>::grid;
        using Diff<TF>::fields;
        using Diff<TF>::boundary;

        enum class Mode {Operational, Diagnostic};
        enum Mom_status : int {
            Mom_invalid_velocity, Mom_zero_velocity, Mom_invalid_shear,
            Mom_invalid_buoyancy, Mom_diagnostic_storage_overflow,
            Mom_invalid_sum, Mom_clipped, Mom_storage_overflow,
            Mom_positive, Mom_capped, Mom_status_count
        };
        enum Scalar_status : int {
            Scalar_zero_gradient, Scalar_invalid_gradient, Scalar_invalid_contraction,
            Scalar_clipped, Scalar_storage_overflow, Scalar_positive,
            Scalar_capped, Scalar_status_count
        };

        struct Aggregates {
            std::uint64_t cells_evaluated = 0;
            std::uint64_t mom[Mom_status_count] {};
            std::uint64_t zero_buoyancy_gradient = 0;
            std::map<std::string, std::vector<std::uint64_t>> scalar;
            double max_evisc = 0.;
            std::map<std::string, double> max_scalar;
        };

        Boundary_cyclic<TF> boundary_cyclic;
        Input& input;
        double dnmax;
        double camd[3];
        double amd_max;
        bool swamd_buoy;
        bool swamd_scalar;
        const Thermo<TF>* thermo = nullptr;
        std::map<std::string, std::string> scalar_coeff;
        std::map<std::string, std::string> scalar_token;
        Aggregates active;
        Aggregates pending;
        Aggregates cumulative;
        bool pending_valid = false;

        void calculate_coefficients(Thermo<TF>&, Mode);
        void fill_coefficient_halos(Field3d<TF>&);
        double timestep_multiplier() const;
        static std::string base32_token(const std::string&);
        static void add_aggregates(Aggregates&, const Aggregates&);
};

#endif
