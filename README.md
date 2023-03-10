
# Style guide
https://www.doc.ic.ac.uk/lab/cplus/cstyle.html

# Decomp & MPI partition

## Decisions

Arbitrary shared-memory comm (`sm_comm`) split where each `sm_comm` handles its own `Send/Recv`.

## Options

1. Arbitrary `sm_comm` split.

    ```
    +-----+-----+-----+    '+' denotes boundary of single proc
    |  1        |  2  |     +--> p_col
    +-----+-----+-----+     |
    |  2  |  3        |     V
    +-----+-----+-----+     p_row
    ```

    1. Leader in each `sm_comm` manages transfers to/from other `sm_comm`.

        Since each `sm_comm` might not have the same relative domain in the cartesian grid, transfers will get very messy.

        In above example, during `Send/Recv` to the right, `sm_comm_1` has to compute the relevant subarrays of `sm_comm_2` and `sm_comm_3`.

    2. Each proc of `sm_comm` manages its own transfer by determining whether its neighbour is internal or external to `sm_comm`.

        In above example, during `Send/Recv` to the right, first proc in `sm_comm_1` only sees one proc in `sm_comm_2`, while second proc in `sm_comm_1` only sees one proc in `sm_comm_3`.

    ## Question
    Are there benefits/drawbacks for splitting `sm_comm_2` into 2 sub-communicators that don't crossover `p_col`? **YES**

    - 2D/3D array indexing will be extremely difficult to handle


2. Enforce size of `sm_comm` to either a) divide `p_col` or b) be divisible by `p_col`.
    ```
    a)                     b)
    +-----+-----+          +-----+-----+
    |  1        |          |  1  |  2  |
    +-----+-----+          +-----+-----+
    |  2        |          |  3  |  4  |
    +-----+-----+          +-----+-----+
    ```
    This can be achieved by:
    - Ensuring `p_col` is set accordingly with `nCPU per node` on particular device to satisfy above conditions.
        - Will need to make sure entire nodes are occupied (i.e. don't share a node with other jobs)
    - Ensuring procs on a node are allocated to the same `p_col`. Options:
        1. Block-allocate the processes, i.e. `(0-23)` on `node 0`, `(24-47)` on `node 1`, ...

            The proc ranks (in `MPI_COMM_WORLD`) on a node **SHOULD** by default ensure `MPI_Cart_create(...,reorder=1,...)` aligns with the shared memory communicator.

        2. Deal-allocate the processes, i.e. `(0,n,2n,...)` on `node 0`, `(1,n+1,2n+1,...)` on `node 1`, ...

            Same as above? Or split communicator by shared memory, then create adjacency matrix for `MPI_Graph_create` to ensure they are aligned in `p_col`.

