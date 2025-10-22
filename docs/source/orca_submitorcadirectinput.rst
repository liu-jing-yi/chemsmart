Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

########################
 Submit Other ORCA Jobs
########################

ChemSmart provides additional ORCA job types for direct input file execution.

**********************************
 Direct ORCA Input File Execution
**********************************

This option allows you to run pre-prepared ORCA input files directly without any modifications by ChemSmart.
Additional files required for the calculation can be specified in the ``inp`` file, and  should be included
in the same directory as the ``.inp`` file.

.. code:: console

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] inp

Basic Usage
===========

-  **Direct execution of ORCA input file**:

      .. code:: console

         chemsmart sub -s shared orca -p direct_run -f my_calculation.inp inp
