module;

#include "cado.h"

#include <cstdint>

#include <array>
#include <string>
#include <type_traits>
#include <vector>
#include <map>

#include "nlohmann/json.hpp"

#include "params.hpp"
#include "timing.h"
#include "verbose.hpp"
#include "fstream_maybe_compressed.hpp"

export module chronograms;

export namespace chronograms /* {{{ */ {
    void configure_switches(cxx_param_list &);
    void declare_usage(cxx_param_list & pl);
    void interpret_parameters(cxx_param_list & pl);
    bool is_enabled();

    // Empty tag structs for parameterless events
    struct NONE {};
    struct INIT {};
    struct QLATTICE {};
    struct SLICING {};
    struct ALLOC {};
    struct AB {};
    struct ECM {};
    struct DUPCHECK {};
    struct BOTCHED {};

    struct SSS   { int side; int level; };
    struct FIB   { int side; int level; uint32_t B; size_t slice; };
    struct DS    { int side; int level; uint32_t B; };
    struct PCLAT { int side; int level; size_t slice; };
    struct PBR   { int M; size_t B; };

    // NOLINTBEGIN(cppcoreguidelines-pro-type-union-access)
    struct bubble_info {
        enum class kind_t : uint8_t {
            NONE,
            INIT,
            QLATTICE,
            SLICING,
            ALLOC,
            SSS,
            FIB,
            AB,
            DS,
            PCLAT,
            PBR,
            ECM,
            DUPCHECK,
            BOTCHED
        };

        static constexpr std::array<const char*, 14> kind_names = {
            "NONE",
            "INIT",
            "QLATTICE",
            "SLICING",
            "ALLOC",
            "SSS",
            "FIB",
            "AB",
            "DS",
            "PCLAT",
            "PBR",
            "ECM",
            "DUPCHECK",
            "BOTCHED", };

        kind_t kind;

        union payload {
            SSS   sss;
            FIB   fib;
            DS    ds;
            PCLAT pclat;
            PBR   pbr;
            DUPCHECK   dupcheck;

            payload() = default;
            explicit payload(SSS e)   : sss(e) {}
            explicit payload(FIB e)   : fib(e) { }
            explicit payload(DS e)    : ds(e) { }
            explicit payload(PCLAT e) : pclat(e) { }
            explicit payload(PBR e)   : pbr(e) { }
        } data = {};

        // NOLINTBEGIN(google-explicit-constructor,hicpp-explicit-conversions)
        bubble_info(NONE)       : kind(kind_t::NONE) {}
        bubble_info(INIT)       : kind(kind_t::INIT) {}
        bubble_info(SSS e)      : kind(kind_t::SSS), data(e) { }
        bubble_info(FIB e)      : kind(kind_t::FIB), data(e) { }
        bubble_info(DS e)       : kind(kind_t::DS), data(e) { }
        bubble_info(PCLAT e)    : kind(kind_t::PCLAT), data(e) { }
        bubble_info(PBR e)      : kind(kind_t::PBR), data(e) { }
	bubble_info(QLATTICE)   : kind(kind_t::QLATTICE) {}
	bubble_info(SLICING)    : kind(kind_t::SLICING) {}
	bubble_info(ALLOC)      : kind(kind_t::ALLOC) {}
	bubble_info(AB)         : kind(kind_t::AB) {}
	bubble_info(ECM)        : kind(kind_t::ECM) {}
	bubble_info(DUPCHECK)   : kind(kind_t::DUPCHECK) {}
	bubble_info(BOTCHED)    : kind(kind_t::BOTCHED) {}
        // NOLINTEND(google-explicit-constructor,hicpp-explicit-conversions)
    };
    // NOLINTEND(cppcoreguidelines-pro-type-union-access)

    struct bubble {
        uint64_t t0 = 0;
        uint64_t t1 = 0;
        uint64_t on_cpu = 0;
        bubble_info info = NONE {};

        static_assert(std::is_trivially_copyable_v<bubble_info>);

        bubble() = default;
        bubble(bubble_info info)
            : info(info)
        {}

        void start() {
            t0 = wct_nanoseconds();
            on_cpu = - microseconds_thread() * 1000;
        }
        void stop() {
            t1 = wct_nanoseconds();
            on_cpu += microseconds_thread() * 1000;
        }
    };

    class [[nodiscard]] bubble_guard {
        std::vector<chronograms::bubble> * destination = nullptr;
        bubble bubble_;

        public:
        static_assert(std::is_trivially_copyable_v<bubble_info>);
        bubble_guard(std::vector<chronograms::bubble> & destination, bubble_info info)
            : destination(&destination)
            , bubble_(info)
        {
            if (!chronograms::is_enabled()) return;
            bubble_.start();
        }

        ~bubble_guard() {
            if (!chronograms::is_enabled() || !destination) return;
            bubble_.stop();
            static_assert(std::is_trivially_copyable_v<chronograms::bubble_info>);
            destination->push_back(bubble_);
        }

        bubble_guard(const bubble_guard&) = delete;
        bubble_guard(bubble_guard&&) = delete;
        bubble_guard& operator=(const bubble_guard&) = delete;
        bubble_guard& operator=(bubble_guard&&) = delete;
        private:
        /* we intentionally make the default ctor private */
        bubble_guard() = default;
        public:
        static bubble_guard dummy() { return {}; }
    };

    std::string format_as(const bubble_info& info);
    void display(std::map<size_t, std::vector<chronograms::bubble>> const & M);

    std::string format_as(const bubble_info& info)
    {
        using kind_t = bubble_info::kind_t;
        auto const & D = info.data;
        // NOLINTBEGIN(cppcoreguidelines-pro-type-union-access)
        switch (info.kind) {
            case kind_t::SSS:
                return fmt::format(
                        "SSS side {} level {}",
                        D.sss.side, D.sss.level);
            case kind_t::FIB:
                return fmt::format(
                        "FIB side {} level {} B {} slice {}",
                        D.fib.side, D.fib.level, D.fib.B, D.fib.slice);
            case kind_t::DS:
                return fmt::format(
                        "DS side {} level {} B {}",
                        D.ds.side, D.ds.level, D.ds.B);
            case kind_t::PCLAT:
                return fmt::format(
                        "PCLAT side {} level {} slice {}",
                        D.fib.side, D.fib.level, D.fib.slice);
            case kind_t::PBR:
                return fmt::format(
                        "PBR M {} B {}",
                        D.pbr.M, D.pbr.B);
            case kind_t::INIT: return "INIT";
            case kind_t::QLATTICE: return "QLATTICE";
            case kind_t::SLICING: return "SLICING";
            case kind_t::ALLOC: return "ALLOC";
            case kind_t::AB: return "AB";
            case kind_t::ECM: return "ECM";
            case kind_t::DUPCHECK: return "DUPCHECK";
            case kind_t::BOTCHED: return "BOTCHED";

            /* we should never see NONE */
            case kind_t::NONE: return "NONE";
        }
        // NOLINTEND(cppcoreguidelines-pro-type-union-access)
        return {};
    }

    void display(std::map<size_t, std::vector<chronograms::bubble>> const & M);
} /* namespace chronograms */ /* }}} */

namespace chronograms {
    static int enable = 0;
    static std::string chronogram_file;

    void configure_switches(cxx_param_list &) {
    }
    void declare_usage(cxx_param_list & pl) {
        pl.declare_usage("chronogram", "Output data to generate a time chart and save it to the given file");
    }
    void interpret_parameters(cxx_param_list & pl) {
        enable = pl.parse("-chronogram", chronogram_file);
    }
    bool is_enabled() { return enable; }

    void display(std::map<size_t, std::vector<chronograms::bubble>> const & M)
    {
        if (!chronograms::is_enabled())
            return;

        size_t total_entries = 0;
        for(auto const & [ k, v ] : M)
            total_entries += v.size();

        verbose_fmt_print (0, 0,
                "# Chronogram info ({} threads, {} entries) will go to {}\n",
                M.size(), total_entries,
                chronograms::chronogram_file);

        ofstream_maybe_compressed out(chronograms::chronogram_file);;

        auto time_min = std::numeric_limits<uint64_t>::max();
        auto time_max = std::numeric_limits<uint64_t>::min();
        auto thr_min = std::numeric_limits<size_t>::max();
        auto thr_max = std::numeric_limits<size_t>::min();

        for(auto const & [ k, v ] : M) {
            for (const auto& b : v) {
                time_min = std::min(time_min, b.t0);
                time_max = std::max(time_max, b.t1);
            }
            thr_min = std::min(thr_min, k);
            thr_max = std::max(thr_max, k);
        }

        verbose_fmt_print (0, 0,
                "# Chronogram info has data for {:.2f} seconds\n",
                double_ratio(time_max - time_min, 1.0e9));

        using json = nlohmann::json;

        json J;
        J["format"] = 20260727;
        J["win_start"] = time_min;
        J["win_end"] = time_max;
        J["categories"] = chronograms::bubble_info::kind_names;
        J["events"] = json();

        for(auto const & [ k, v ] : M) {
            for (const auto& b : v) {
                json E;
                E.push_back(b.t0 - time_min);
                E.push_back(b.t1 - b.t0);
                E.push_back(static_cast<int>(b.info.kind));
                E.push_back(b.on_cpu);
                E.push_back(fmt::format("{}", b.info));
                J["events"][std::to_string(k)].push_back(E);
            }
        }
        out << J;

        verbose_fmt_print (0, 0,
                "# Chronogram info can be viewed"
                " using https://cado-nfs.inria.fr/chronogram.html"
                " or ./scripts/chronograms/chronograms.html\n");
    }
}
