#include "cado.h" // IWYU pragma: keep
#include <ostream>    // std::ostream // IWYU pragma: keep
#include <iomanip> // for std::fixed and std::setprecision
#include <pthread.h>
#include <sys/time.h>
#include "verbose.h"
#include "tdict.hpp"
#include "params.h"

int time_bubble_chaser::enable = 0;

namespace tdict {

/* default value is set in configure_switches */ 
int global_enable;

#ifndef DISABLE_TIMINGS

pthread_mutex_t slot_base::m = PTHREAD_MUTEX_INITIALIZER;

void declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "T",   "enable fine-grain timings (use twice to get them for each q)");
}

void configure_aliases(cxx_param_list &)
{
}

void configure_switches(cxx_param_list & pl)
{
    param_list_configure_switch(pl, "-T", &global_enable);
    /* We now rely on fine-grain switches to provide _all_ timings */
    global_enable = 1;
}

std::ostream& operator<<(std::ostream & o, timer_seconds_thread_and_wct::type const & a) {
    o << a.t << " (" << a.w << " wct";
    if (a.w) o << ", " << a.t/a.w*100.0 << "% cpu";
    o << ")";
    return o;
}

#else

void declare_usage(cxx_param_list &) {}
void configure_aliases(cxx_param_list &) {}
void configure_switches(cxx_param_list &) {}

#endif

};

time_bubble_chaser::time_bubble_chaser(int thread, kind_t kind, id_t const & id): thread(thread), kind(kind), id(id) {
    if (!enable) return;
    gettimeofday(&tv_get, NULL);
    on_cpu = - (double) microseconds_thread();
}
time_bubble_chaser& time_bubble_chaser::put() {
    if (!enable) return *this;
    gettimeofday(&tv_put, NULL);
    on_cpu += microseconds_thread();
    return *this;
}
void timetree_t::display_chart() const
{
    verbose_output_print (0, 2, "# time chart has %zu entries\n", chart.size());
    for(auto const & T : chart) {
        std::ostringstream os;

        os << "thread " << T.thread;

        /* This switch should really be passed as a function pointer */
        switch(T.kind) {
            case time_bubble_chaser::FIB:
                os << " FIB"
                    << " side " << T.id[0]
                    << " level " << T.id[1]
                    << " B " << T.id[2]
                    << " slice " << T.id[3];
                break;
            case time_bubble_chaser::DS:
                os << " DS"
                    << " side " << T.id[0]
                    << " level " << T.id[1]
                    << " B " << T.id[2];
                break;
            case time_bubble_chaser::SSS:
                os << " SSS"
                    << " side " << T.id[0]
                    << " level " << T.id[1];
                break;
            case time_bubble_chaser::AB:
                os << " AB";
                for(auto const & x : T.id)
                    os << " x" << (&x-begin(T.id)) << " " << x;
                break;
            case time_bubble_chaser::PBR:
                os << " PBR";
                os  << " M " << T.id[1]
                    << " B " << T.id[2];
                //for(auto const & x : T.id)
                    //os << " x" << (&x-begin(T.id)) << " " << x;
                break;
            case time_bubble_chaser::PCLAT:
                os << " PCLAT"
                    << " side " << T.id[0]
                    << " level " << T.id[1]
                    << " slice " << T.id[3];
                break;
            case time_bubble_chaser::ECM: os << " ECM"; break;
            case time_bubble_chaser::INIT: os << " INIT"; break;
            case time_bubble_chaser::BOTCHED: os << " BOTCHED"; break;
            case time_bubble_chaser::SKEWGAUSS: os << " SKEWGAUSS"; break;
            case time_bubble_chaser::SLICING: os << " SLICING"; break;
            case time_bubble_chaser::ALLOC: os << " ALLOC"; break;
            case time_bubble_chaser::ADJUST: os << " ADJUST"; break;
        };

        double t0 = T.tv_get.tv_sec + 1.0e-6 * T.tv_get.tv_usec;
        double t1 = T.tv_put.tv_sec + 1.0e-6 * T.tv_put.tv_usec;
        double t = T.on_cpu;
        os  << std::fixed << std::setprecision(9)
            << " t0 " << t0
            << " t1 " << t1
            << " time " << t;

        verbose_output_print (0, 2, "#  %s\n", os.str().c_str());
    }
}


template class std::map<tdict::key, tdict::slot_base const *>;

template struct tdict::tree<tdict::timer_seconds_thread>;
template class std::map<tdict::key, tdict::tree<tdict::timer_seconds_thread> >;
// template struct std::pair<tdict::key const, tdict::slot_base const *>;
template struct tdict::tree<tdict::timer_seconds_thread>::accounting_child_meta<tdict::tree<tdict::timer_seconds_thread>::accounting_base>;

#ifdef  HAVE_GCC_STYLE_AMD64_INLINE_ASM
template struct tdict::tree<tdict::timer_ticks>;
template class std::map<tdict::key, tdict::tree<tdict::timer_ticks> >;
template struct tdict::tree<tdict::timer_ticks>::accounting_child_meta<tdict::tree<tdict::timer_ticks>::accounting_base>;
#else
template struct tdict::tree<tdict::timer_none>;
template class std::map<tdict::key, tdict::tree<tdict::timer_none> >;
template struct tdict::tree<tdict::timer_none>::accounting_child_meta<tdict::tree<tdict::timer_none>::accounting_base>;
#endif

