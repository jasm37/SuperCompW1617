# How to run Intel Analyzer
  1.  First make sure that you are tunneling through lxhalle from you local to make the X11 display work (for the tracer UI as Totalview).
  2.  Unload the module **mpi.ibm** and load **mpi.intel** and **itac**.
  3.  You do not need to recompile the code, just add the **-trace** in the mpiexec command in the Load Leveler script:

      ```sh
      mpirun -trace -n 32 ./gauss ./ge_data/size512x512
      ```

      This will produce many stf files in the **/ge_data** folder named **gauss.stf..**. You visualize them by running:

      ```sh
      traceanalyzer gauss.stf
      ```

  Source: https://www.lrz.de/services/software/parallel/itac/#usage
