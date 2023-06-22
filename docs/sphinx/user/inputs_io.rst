.. _inputs_io:

Section: io
~~~~~~~~~~~~~~~~~

This section deals with parameters that control input/output to the simulation.

.. input_param:: io.check_file

   **type:** String, optional, default = "chk"

   If :input_param:`time.checkpoint_interval` is greater than zero this is the name of the checkpoint 
   file appended with the current timestep
   
.. input_param:: io.plot_file

   **type:** String, optional, default = "plt"

   If :input_param:`time.plot_interval` is greater than zero this is the name of the plot
   file appended with the current timestep
   
.. input_param:: io.restart_file

   **type:** String, optional, default = ""

   If a string is present `amr-wind` will restart using the specified file in the string.
   
   
