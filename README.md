Please not on some server, we need to change the compiler link command at the end of lines
1. The program takes the features of MPI_IO. the program read the whole trajectories first then process. This enables faster processing, but will consumes a lot of memory.
2. It could also read one snapshot and process, then the other one. This will result in slower processing but smaller memeories.
3. The program is written and developed by **Jiahao Zhang** under tne supervision of **Prof.Rappe**
