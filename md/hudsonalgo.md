---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Hudson's linear-time simulation algorithm.

In a 1990, Dick Hudson published the C code for a linear time algorithm to simulate gene genealogies under the Kingman coalescent {cite}`Hudson1990-ff`.
His method works by realizing that for a sample of $n$ chromosomes, the final tree with have $2n - 1$ nodes that we may label $[0, 2n)$. 
The first $n$ nodes correspond to the present day sample, and the remaining $n-1$ are the possible ancestors.

To get a valid tree, we require:

* The value at the $i^{th}$ element of the array is the index of the node ancestral to $i$.
* The ancestor of the root node is a `NULL` value, for which we will use $-1$.

For example, if the common ancestor of sample nodes 1 and 3 is node 11, then `tree[1]` and `tree[3]` will both store the value 11.
The root node is the final value in the array, and `tree[root]` will be -1, signalling "this node has no parent".

More generally, at any time in the past where $i$ lineages remain,
the first $i$ nodes are represent the sample, and the remaining nodes
are their possible parents. We wish to take a random pair from the first $i$
nodes and assign as their parent the first of the remaining possible ancestral nodes.

The algorithm is:

A. *Initialization:* for a sample of size $n$:
   * Initialize a tree, $T$, to $2n - 1$ `NULL` values of $-1$.
   * Initialize a list of times, $T_n$, to $2n - 1$ `NULL` values of "not a number".
   * Initialize a list of remaining nodes, $R$, to contain all values $[0, n)$.
   * Set $i = n$, $t = 0$, and $X = 2\times n$.

B. Termination:
   * If $i = 1$, end the algorithm
   * If $i \gt 1$, continue to step `C`.
   
C. *Determine the ancestor*:
   * Set ancestor, $A$ equal to $X - i$.

D. *Determine coalescence time*:
   * Generate an exponential random deviate, $e$, using rate ${i \choose 2},$ which equals $\frac{i(i-1)}{2}.$

E. *Update the current time*:
   * Set $t$ equal to $t + e$.

F. *Update the list of times*:
   * Set $T_c[A]$ equal to $t$.
    
G. *Coalescence*: 
   * Pick two values, $j$ and $k$, uniformly and without replacement from $[0, i)$.
   * Set $c_1$ and $c_2$ equal to the $j^{th}$ and $k^{th}$ values in $R$, respectively.
   * Set $T[c_1]$ and $T[c_2]$ equal to $A$.
   * Set $R[min(j,k)]$ to $A$
   * Set $R[max(j,k)]$ to $R[i-1]$
   * Set $i$ to $i - 1$

H. *Continue*:
   * Return to step `B`.

The updating of $R$ in step `G` is the trick.

The above algorithm is very straightforward to write out in Python.
We use `numpy` arrays for efficiency.

```{code-cell}
import typing
import numpy as np

def h1990(nsam: int) -> typing.Tuple[np.ndarray, np.ndarray]:
    """
    The linear-time algorithm of Hudson, 1990.

    The citation for this algorithm is
    Hudson, Richard R. 1990.
    “Gene Genealogies and the Coalescent Process.”
    Oxford Surveys in Evolutionary Biology 7 (1): 44.

    Time is scaled in units of 2N generations.

    :param nsam: The sample size
    :type nsam: int
    
    :return: The tree and the times
    :rtype: tuple
    """
    tree = np.array([-1]*(2*nsam-1), dtype=np.int32)
    times = np.array([np.nan]*len(tree))
    times[0:nsam][:] = np.zeros(nsam)
    
    time = 0.0
    i = nsam
    R = np.arange(n, dtype=np.int32)
    while i > 1:
        # Generate time to next coalescent event,
        # in units of 2N generations.
        rcoal = (i*(i-1))/2.
        tcoal = np.random.exponential(1./rcoal)
        time += tcoal

        # This is the index of the
        # ancestor node
        A = 2*nsam - i
        times[A] = time

        # Perform the swap steps
        # of the algorithm
        c = np.random.choice(i, 2, replace=False)
        c1 = c[0]
        c2 = c[1]
        n1 = R[c1]
        n2 = R[c2]
        tree[n1] = A
        tree[n2] = A
        R[min(c1,c2)] = A
        R[max(c1,c2)] = R[i-1]
        
        i -= 1

    return (tree, times)
```

## Applications

```{code-cell}
n = 10
nreps = 5000
np.random.seed(54321)
```

### Time to the most recent common ancestor

```{code-cell}
tmrca = 0.0

for _ in range(nreps):
    tree, times = h1990(n)
    # What am I doing in the following line? 
    tmrca += times[-1]
print(tmrca/nreps)
```

The expectation is $2\left(1-\frac{1}{n}\right)$, which you can get from Wakeley, chapter 3, equation 3.24 and lots of other places:

```{code-cell}
2*(1 - 1/n)
```

### Total time on the tree

Given the tree structure, we can get the total time using a linear algorithm.

```{code-cell}
ttimes = np.zeros(nreps)

def get_ttime(tree: np.ndarray, times: np.ndarray) -> float:
    ttime = 0.0
    for i, node in enumerate(tree):
        if node != -1:
            branch_len = times[node] - times[i]
            assert branch_len > 0.0, f"{branch_len} {i} {node} {tree}"
            ttime += branch_len
    return ttime
            

for i in range(nreps):
    tree, times = h1990(n)
    ttimes[i] = get_ttime(tree, times)
```

We get the following mean and variance from the simulations:

```{code-cell}
print(f"The mean total time is {ttimes.mean():.4f}")
print(f"The variance in total time is {ttimes.var():.4f}")
```

The analytical expectation and variance are $4\sum_{i=1}^{n-1}\frac{1}{i}$ and $4\sum_{i=1}^{n-1}\frac{1}{i^2},$ respectively:

```{code-cell}
def ttime_theory(n: int) -> typing.Tuple[float, float]:
    ett = 0.
    for i in range(1, n):
        ett += 1/i
    ett *= 2.0
    vtt = 0.
    for i in range(1, n):
        vtt += 1/(i**2)
    vtt *= 4.0
    return ett, vtt
    
print(f"The mean and variance for n={n} are {ttime_theory(n)}")
```

### The expectation of the "site frequency spectrum"

The expected amount of time leading to $i$ tree tips is $2/i$, in units of $2N_e$ generations.
We can get these times using an efficient tree traversal from tips to root:

```{code-cell}
def update_sfs_times(nsam: int, tree: np.ndarray, times: np.ndarray, sfs_times: np.ndarray) -> np.ndarray:
    ndescendants = np.zeros(len(tree), dtype=np.int32)
    for i in range(nsam):
        p = i 
        while p != -1:
            ndescendants[p] += 1
            p = tree[p]
    
    for i, n in enumerate(ndescendants):
        if tree[i] != -1:
            sfs_times[n-1] += times[tree[i]] - times[i]
    return sfs_times
        

sfs_times = np.zeros(n - 1)
for i in range(nreps):
    tree, times = h1990(n)
    sfs_times = update_sfs_times(n, tree, times, sfs_times)

sfs_times /= nreps
sfs_times
```

The expectation for these bins is:

```{code-cell}
2./(np.arange(n-1)+1.0)
```

Now, we have the expected *time* for each bin of the `SFS`.
If $\theta = 4N_e\mu$ is our scaled mutation rate, then $\theta/2$ times each entry in the above array gives the expected number of mutations in each bin.
