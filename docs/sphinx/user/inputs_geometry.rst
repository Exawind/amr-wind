Section: geometry
~~~~~~~~~~~~~~~~~~~~~

This section deals with inputs related to the problem domain.

.. input_param:: geometry.prob_lo

   **type:** List of 3 real numbers, mandatory

   The coordinates of *lower corner* of the computational domain bounding box.

.. input_param:: geometry.prob_hi

   **type:** List of 3 real numbers, mandatory

   The coordinates of the *upper corner* of the computational domain bounding box.

.. input_param:: geometry.is_periodic

   **type:** List of 3 integers, mandatory

   Flags indicating whether the flow is periodic in the ``x``, ``y``, or ``z``
   directions respectively.
