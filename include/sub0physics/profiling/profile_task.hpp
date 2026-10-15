#pragma once

#ifndef SUB0PHYSICS_ITT
#define SUB0PHYSICS_ITT 0
#endif

#if SUB0PHYSICS_ITT
#include <ittnotify.h>
#endif

namespace sub0physics
{

/**
 * Marks a named task on the calling thread for Intel VTune and other ITT-aware profilers.
 *
 * With SUB0PHYSICS_ITT off (the default) it compiles to nothing. With it on, VTune's timeline and
 * "Task" grouping show where each run, phase or batch starts and ends. ITT interns domain and task
 * names, so creating one is a lookup, not an allocation; still, mark coarse phases such as a whole run
 * or a recording interval, never a single integration step.
 *
 * @note Not copyable or movable; begins in the constructor and ends in the destructor on one thread.
 */
class ProfileTask
{
public:
    /**
     * Begins a task.
     *
     * @param domain Grouping name, e.g. the library; string literal.
     * @param name Task name; string literal.
     */
    ProfileTask([[maybe_unused]] const char* domain, [[maybe_unused]] const char* name) noexcept
    {
#if SUB0PHYSICS_ITT
        domain_ = __itt_domain_create(domain);
        __itt_task_begin(domain_, __itt_null, __itt_null, __itt_string_handle_create(name));
#endif
    }

    ~ProfileTask()
    {
#if SUB0PHYSICS_ITT
        __itt_task_end(domain_);
#endif
    }

    ProfileTask(const ProfileTask&) = delete;
    ProfileTask& operator=(const ProfileTask&) = delete;
    ProfileTask(ProfileTask&&) = delete;
    ProfileTask& operator=(ProfileTask&&) = delete;

private:
#if SUB0PHYSICS_ITT
    __itt_domain* domain_; // non-owning; ITT owns domains for the process lifetime
#endif
};

} // namespace sub0physics
