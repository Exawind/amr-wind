.. _inputs_ransconvergence:

Section: RANSConvergence
~~~~~~~~~~~~~~~~~~~~~~~~

This section controls the pseudo-steady convergence monitor. The monitor
samples horizontal wind speed and turbulent kinetic energy at a list of
user-supplied points, and once every point has stopped changing it asks the
time integrator to stop, which writes a final plotfile and checkpoint. It
removes the need to guess a :input_param:`time.stop_time` that is long enough
to reach a steady state but not wastefully longer.

The prefix is the label set in ``incflo.post_processing``. For example
``incflo.post_processing = convergence``. The inputs controlling output timing
that are shared with other post-processing types are listed in the
[:ref:`post-processing section <inputs_post_processing>`]; this monitor drives
its own sampling cadence and does not use them.

Restricted to the KLAxell model
```````````````````````````````

The monitor aborts unless :input_param:`turbulence.model` is ``KLAxell``.

The criterion is a drift test on an instantaneous sampled value, which is only
meaningful when that value approaches a limit. A RANS run marching to a
pseudo-steady state has such a limit. A large eddy simulation does not: the
velocity at a point stays turbulent indefinitely, so a drift test on the raw
value never passes, and a drift test on a running mean passes trivially because
a cumulative mean converges as :math:`1/N` whether or not the flow is
stationary. Applying the monitor to an LES would report a convergence the
solution does not have, so it is refused rather than offered with a warning.

How convergence is measured
```````````````````````````

At each sample the monitor interpolates velocity and ``tke`` at every monitor
point and stores them. Over a trailing window of width
:input_param:`convergence.window` it then measures the peak-to-trough spread of
horizontal wind speed :math:`\sqrt{u^2+v^2}` and of ``tke``. A point has
converged when both spreads are below tolerance, and the run stops when every
point has converged.

The spread over a window is used rather than a difference between successive
samples because an atmospheric boundary layer under Coriolis and geostrophic
forcing does not decay monotonically to steady state. It spirals in through a
damped inertial oscillation of period :math:`2\pi/f`, which for
:math:`f = 10^{-4}\,\mathrm{s^{-1}}` is about 17.5 hours. A difference between
successive samples vanishes at every turning point of that oscillation, so a
monitor built on it would report convergence at the first peak, hours before
the oscillation has actually decayed. A window still contains the swing on
both sides of a turning point and is not fooled.

A window shrinks the turning-point problem but does not remove it. Near a
turning point the signal is locally quadratic, so the spread measured over a
window of width :math:`W` falls to roughly :math:`W/T` of its typical value,
where :math:`T` is the oscillation period. On a one-hour window with a 17.5
hour inertial period, that is a dip of more than an order of magnitude: in a
test case the speed envelope peaked at 0.52 m/s and dipped to 0.022 m/s at a
turning point, while the flow was still far from steady.

The remedy is :input_param:`convergence.hold_time`. The dip lasts on the order
of one window, so requiring the criterion to hold continuously for longer than
that means a turning point cannot stop the run whatever the window and
tolerance are. The default of two windows is measured against that dip and
should not be lowered without a reason.

Three consequences for choosing inputs:

* :input_param:`convergence.window` should be a meaningful fraction of the
  inertial period, one to two hours for a mid-latitude or polar case. Longer
  windows make the turning-point dip shallower.
* :input_param:`convergence.hold_time` guards the dip that the window leaves
  behind. Leave it at the default unless you have measured the dip in your own
  configuration.
* :input_param:`convergence.sample_interval_time` should be generous. Sampling
  every step buys no information, because the solution cannot move
  meaningfully in one timestep.

The vertical velocity is deliberately excluded from the speed. In a converging
boundary layer it is small and comparatively noisy, so including it adds
jitter to the envelope and delays the stop without making the test more
informative.

Every point is tested independently and all points must pass. The metric is
never averaged across points, because a point drifting up and a point drifting
down would cancel in the mean and the aggregate would pass while neither point
had converged.

Choosing tolerances
```````````````````

Each quantity takes an absolute and a relative tolerance, and the test uses
whichever is larger:

.. math::

   \mathrm{spread} < \max(\mathrm{abs\_tol},\;
                          \mathrm{rel\_tol} \cdot |\overline{\phi}|)

A single absolute tolerance does not work for ``tke``, which spans orders of
magnitude between the surface layer and the region above the boundary-layer
top. A single relative tolerance fails at the other end, becoming
unsatisfiable as ``tke`` approaches zero aloft. The combined form behaves as
an absolute floor where the quantity is small and as a relative test where it
is large.

Both absolute tolerances must be greater than zero and neither relative
tolerance may be negative; the monitor refuses to start otherwise.

Start by running with :input_param:`convergence.stop_on_convergence` set to
false on a case whose convergence you already trust. The monitor then reports
the spreads at each check without ever stopping the run, which is the cheapest
way to see the noise floor of your configuration and pick tolerances that are
neither unreachable nor met by accident.

:input_param:`time.max_step` and :input_param:`time.stop_time` remain in force
as backstops. A tolerance set below the noise floor of the run will simply
never trigger.

Estimated time to convergence
`````````````````````````````

When the envelope is decaying roughly exponentially, the monitor fits
:math:`s(t) = A e^{-rt}` by least squares to the recent history of the worst
normalized spread and extrapolates to the point where it would reach
tolerance. The result is printed as an estimated time to convergence, which is
useful for deciding whether a queued job has enough wall time left to finish.

The estimate is a diagnostic and never a stopping criterion: the run stops
only on measured spreads. Extrapolating an exponential fit to a threshold is
precisely where such a fit is least reliable, so the estimate is suppressed
unless the fit is good, controlled by :input_param:`convergence.eta_min_rsq`.
It is also suppressed when the envelope is not decaying, which is the honest
report when a run has stalled on a noise floor rather than a confident
prediction that will never come true.

Monitor points
``````````````

Points are read from a text file whose first line is the number of points,
followed by one ``x y z`` triple per line. This is the same format used by
:ref:`ProbeSampler <inputs_sampling>`.

The monitor prints the coordinates of every point it is monitoring at startup.
Check that list: a point outside the domain is silently moved inside it by the
underlying probe sampler, with a warning, so the printed coordinates are the
ones actually being used.

If a terrain height field is present, points below the terrain surface are a
fatal error. Such a point is held near zero by the drag forcing and converges
immediately, and because every point must pass the test it would weaken the
criterion without ever failing.

Restarts
````````

The window is not carried across a restart. After a restart the monitor
refills the window from scratch before it can declare convergence, which costs
at most one window of extra runtime. This is deliberate: a restart may change
resolution, forcing or terrain, and carrying stale samples across it would let
a run stop on evidence gathered under different physics.

Inputs
``````

.. input_param:: convergence.type

   **type:** String, mandatory

   To use the convergence monitor specify with keyword ``RANSConvergence``

.. input_param:: convergence.probe_location_file

   **type:** String, optional, default = ``probe_locations.txt``

   Path to the file listing the monitor points.

.. input_param:: convergence.start_time

   **type:** Real, mandatory

   Simulation time at which monitoring begins. There is no default and a value
   of zero is rejected. A KLAxell run is driven toward its target by forcing
   terms that ramp in over time, such as
   :input_param:`DragForcing.bc_forcing_time_factor` and the mesoscale sponge.
   Convergence measured before those ramps have finished says only that the
   solution is tracking a moving target, so this delay must be chosen
   deliberately rather than defaulted.

.. input_param:: convergence.sample_interval_time

   **type:** Real, mandatory

   Simulation time in seconds between samples. Specified in physical time
   rather than timesteps because these runs use adaptive timestepping.

.. input_param:: convergence.window

   **type:** Real, mandatory

   Width in seconds of the trailing window over which the peak-to-trough
   spread is measured. Must be larger than the sample interval.

.. input_param:: convergence.velocity_abs_tol

   **type:** Real, optional, default = 0.01

   Absolute tolerance in m/s on the horizontal wind speed spread. Must be
   greater than zero: it is the floor that keeps the test meaningful where the
   speed is small, every check divides a spread by the resulting tolerance, and
   a value of zero could never be met in any case.

.. input_param:: convergence.velocity_rel_tol

   **type:** Real, optional, default = 0.001

   Relative tolerance on the horizontal wind speed spread, applied to the
   window mean of the speed at that point.

.. input_param:: convergence.tke_abs_tol

   **type:** Real, optional, default = 0.001

   Absolute tolerance in m^2/s^2 on the turbulent kinetic energy spread. Must
   be greater than zero, for the same reason as the velocity floor above.

.. input_param:: convergence.tke_rel_tol

   **type:** Real, optional, default = 0.01

   Relative tolerance on the turbulent kinetic energy spread, applied to the
   window mean of ``tke`` at that point.

.. input_param:: convergence.hold_time

   **type:** Real, optional, default = twice ``window``

   How long in seconds every point must stay within tolerance continuously
   before the run is stopped. A single passing check is not enough, because the
   envelope dips at every turning point of the inertial oscillation. The
   elapsed hold is reported in the output file, and it resets to zero as soon
   as any point falls back outside tolerance.

.. input_param:: convergence.min_samples

   **type:** Integer, optional, default = 4

   Fewest samples that must be present in the window before convergence may be
   declared. Must be at least 2 for a spread to be meaningful.

.. input_param:: convergence.stop_on_convergence

   **type:** Boolean, optional, default = true

   When true, the run stops once every point has converged, and a plotfile and
   checkpoint are written unconditionally. When false the monitor reports at
   each check but never stops the run.

.. input_param:: convergence.report_eta

   **type:** Boolean, optional, default = true

   Report an estimated time to convergence from an exponential fit to the
   envelope history.

.. input_param:: convergence.eta_fit_window

   **type:** Real, optional, default = 10 times ``window``

   Trailing span in seconds of envelope history used for the exponential fit.
   Using the whole history instead would let the early transient, before the
   forcing ramps have settled, contaminate the fitted decay rate.

.. input_param:: convergence.eta_min_samples

   **type:** Integer, optional, default = 5

   Fewest envelope samples that will be fitted.

.. input_param:: convergence.eta_min_rsq

   **type:** Real, optional, default = 0.5

   Smallest coefficient of determination, computed in log space, for which the
   estimated time to convergence is reported. A fit below this is discarded
   rather than quoted.

Output
``````

A text file named after the label is written to the post-processing directory,
with one line per check recording the time, how many samples the window holds
and whether it is full, how many points have converged, the worst point for
each quantity together with its spread and the tolerance it was compared
against, how long the criterion has held, and the estimated time to
convergence. The same information is printed to the log.

Read ``window_full`` before ``num_converged``. While the window is still
filling it holds only a couple of samples, whose spread is small whatever the
flow is doing, so early rows routinely show every point converged without
meaning it.

The worst-offender columns are the ones to read when a run does not converge:
they distinguish a flow that is genuinely still evolving from a single badly
placed point holding up an otherwise converged field, and from a tolerance
that is simply too tight.

Example
```````

::

   incflo.post_processing            = convergence

   convergence.type                  = RANSConvergence
   convergence.probe_location_file   = monitor_points.txt

   convergence.start_time            = 7200.0
   convergence.sample_interval_time  = 60.0
   convergence.window                = 3600.0

   convergence.velocity_abs_tol      = 0.01
   convergence.velocity_rel_tol      = 0.001
   convergence.tke_abs_tol           = 0.001
   convergence.tke_rel_tol           = 0.01

   convergence.stop_on_convergence   = true
