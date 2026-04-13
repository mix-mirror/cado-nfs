#include "cado.h" // IWYU pragma: keep

#include <cerrno>
#include <cstdio>
#include <cstring>

#include <string>
#include <vector>
#if !defined(__x86_64) && !defined(__i386)
#include <atomic>
#endif

#include <pthread.h>
#ifdef HAVE_GETRUSAGE
#include <sys/resource.h>
#endif
#include <gmp.h>

#include "cado_popen.h"
#include "macros.h"
#include "ringbuf.hpp"
#include "portability.h"
#include "timing.h"
#include "filter_io.hpp"

void * cado::filter_io_details::filter_rels_producer_thread(
    ringbuf & r,
    std::vector<std::string> const & input_files,
    timingstats_dict_ptr stats)
{
    for(auto const & filename : input_files) {
        int status;
        /* We expect all the "filenames" to have been returned by
         * prepare_grouped_command_lines, thus in fact be commands to be
         * passed through popen()
         */
        FILE * f = cado_popen(filename.c_str(), "r");
        ssize_t const rc = r.feed_stream(f);
        int const saved_errno = errno;
#ifdef  HAVE_GETRUSAGE
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
        struct rusage rus;
        status = cado_pclose2(f, &rus);
        if (stats) timingstats_dict_add(stats, "feed-in", &rus);
#else
        status = cado_pclose(f);
#endif
        if (rc < 0) errno = saved_errno;

        if (rc < 0 || status == -1
#if defined(WIFEXITED) && defined(WEXITSTATUS)
            || !WIFEXITED(status) || WEXITSTATUS(status) != 0
#endif
           ) {
            fprintf(stderr,
                    "%s: load error (%s) while %s\n%s\n",
                    __func__,
                    strerror(errno),
                    rc < 0 ? "reading from" : "closing",
                    filename.c_str());
            abort();
        }
    }
    r.mark_done();
    if (stats) timingstats_dict_add_mythread(stats, "producer");
    /*
    double thread_times[2];
    thread_seconds_user_sys(thread_times);
    fprintf(stderr, "Producer thread ends after having spent %.2fs+%.2fs on cpu\n",
            thread_times[0],
            thread_times[1]);
            */
    return nullptr;
}
